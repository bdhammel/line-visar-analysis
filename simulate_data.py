import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp2d
from scipy.ndimage import convolve1d
from PIL import Image

c = 299792.0  # um/ns


class Ray:

    def __init__(self, lambda0: "um" = .532,
                 pulse_length: "ns" = 10, radius: "um" = 100):
        """
        Parameters
        ----------
        lambda_: flaot
            wave length of the light in um
        pulse_length : float
            time in us
        radius : float
            radius of the beam
        """
        self.radius = radius
        self.lambda0 = lambda0
        self._t = np.linspace(0, pulse_length, 2048)
        self._y = np.linspace(-2, 2, 2048)
        self._tt, self._yy = np.meshgrid(self._t, self._y)
        self._lambda0 = lambda0*np.ones_like(self._tt)
        self._delta = 0

    @property
    def pulse_length(self):
        return self._t

    @property
    def beam_width(self):
        return self._y

    @property
    def _k(self):
        return 2*np.pi/self._lambda

    @property
    def _k0(self):
        return 2*np.pi/self._lambda0

    @property
    def phi(self):
        return self._k*self.dz + self._k0*self.dz

    def E(self, t):
        E = np.exp(1j*(self.phi + self._delta))
        E_real = np.real(E)
        E_imag = np.imag(E)
        fE_real = interp2d(self._t, self._y, E_real)
        fE_imag = interp2d(self._t, self._y, E_imag)
        return fE_real(t, self._y) + 1j*fE_imag(t, self._y)

    def set_lambda(self, lambda_):
        self._lambda = lambda_

    def propogate(self, dz):
        self.dz = dz

    def add_delta(self, delta):
        self._delta = delta


class Target:

    def __init__(self, velocity_equation):
        """
        Parameters
        ----------
        velocity_equation : str or fn
            either step or sigmoid to use default velocity profile, or a
            function that excepts a t and y meshgrid
        """
        self._t = np.linspace(-5, 15, 2048)
        self._y = np.linspace(-3, 3, 2048)
        self.tau = 0

        self._tt, self._yy = np.meshgrid(self._t, self._y)
        if velocity_equation == "step":
            self.velocity_equation = self.step
        elif velocity_equation == "sigmoid":
            self.velocity_equation = self.sigmoid
        elif velocity_equation == "stationary":
            self.velocity_equation = self.stationary
        else:
            self.velocity_equation = velocity_equation

    @property
    def _zz(self):
        dt = np.diff(self._t).mean()
        return np.cumsum(self._vv, axis=1)*dt

    @property
    def zz(self):
        return interp2d(self._t, self._y, self._zz)

    @property
    def _dz(self):
        """Path the light travels to the target and back
        """
        dzz = self._zz[..., -1, np.newaxis] - self._zz
        return dzz

    @property
    def dz(self):
        return interp2d(self._t, self._y, self._dz)

    @property
    def _vv(self):
        return self.velocity_equation(self._tt, self._yy)

    @property
    def vv(self):
        return interp2d(self._t, self._y, self._vv)

    @staticmethod
    def sigmoid(t: "ns", y: "um", max_velocity: "um/ns" = 5):
        """A velocity profile that follows a sigmoid like shape
        """
        return max_velocity*np.exp(-y**4)/(1 + np.exp(-5*(t-3)))

    @staticmethod
    def step(t: "ns", y: "um", max_velocity: "um/ns" = 1):
        """A discontinuous jump velocity profile
        """
        assert t.shape == y.shape
        v = np.zeros_like(t)
        v[t > 3] = max_velocity
        return v

    @staticmethod
    def stationary(t: "ns", y: "um"):
        """A static target, not moving
        """
        return np.zeros_like(t)

    def reflect_off_target(self, ray):
        ray = self._doppler_shift(ray)
        dz = self.dz(ray.pulse_length, ray.beam_width)
        ray.propogate(dz)
        return ray

    def _doppler_shift(self, ray):
        vv = self.vv(ray.pulse_length, ray.beam_width)
        ray.set_lambda(ray.lambda0*(1 - 2*vv/c))
        return ray

    def reflection_intensity(self, ray):
        dy = np.diff(ray.beam_width).mean()
        dz = np.diff(self.zz(ray.pulse_length, ray.beam_width), axis=0)
        theta = np.arctan(dz/dy)
        Idot = np.vstack(
            (np.ones(shape=(2048)), np.apply_along_axis(np.cos, 0, theta))
        )
        return Idot

    def plot_velocity(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        Axes3D.plot_surface(ax, self._tt, self._yy, self._vv)
        ax.set_xlabel("Time [ns]")
        ax.set_ylabel("x [mm]")
        ax.set_zlabel("Velocity [km s-1]")

        fig = plt.figure()
        im = plt.pcolormesh(self._tt, self._yy, self._vv)
        cb = fig.colorbar(im)
        plt.xlabel("Time [ns]")
        plt.ylabel("x [mm]")
        cb.set_label("Velocity [km s-1]")


class Etalon:

    def __init__(self, thickness: "mm", n):
        """Initial an etalon object

        Parameters
        -----------
        d : float
            thickness of the etalon
        n : float
            ndex of refraction of the etalon
        """
        self._n = n
        self._d = thickness

    @property
    def tau(self) -> "ns":
        return 2*self._d/c*(self._n - 1/self._n)

    def VPF(self, lambda0=.532):
        return lambda0/(2*self.tau)

    def set_VPF(self, VPF, lambda0: "um"):
        tau = lambda0/(2*VPF)
        self.set_tau(tau)

    def set_tau(self, tau: "ns"):
        """Change the thickness of the etalon to match a
        """
        self._d = c*tau/(2*(self._n - 1/self._n))


class Interferometer:

    def __init__(self, etalon, tau: "ns" = .1):
        """

        Parameters
        ----------
        etalon : Etalon
            the etalon used in the interferometer, provides VPF
        tau : float
            the time resolution of the streak camera, determined by the
            width of the streak slit
        """
        self.etalon = etalon
        self.tau = tau

    def _interfear_ray(self, ray):
        """Generate the interference pattern
        """
        # get the electric field over the pulse length
        E1 = ray.E(ray.pulse_length)
        # generate the offset for the second ray
        _delta_shape = len(ray.beam_width)
        ray.add_delta(
            np.linspace(0, 100, _delta_shape).reshape(_delta_shape, 1)
        )

        # generate the second ray, which is delayed by the etalon thickness
        E2 = ray.E(ray.pulse_length - self.etalon.tau)

        # Super position of the rays
        E = E1 + E2

        # only take the real component of the inner product (intensity)
        Icos = np.real(E*E.conj())

        return Icos

    def _add_noise(self, im, ray, target, noise_level, signal_level):
        """Add detector noise to the generated fringe pattern
        """

        print("...Including noise")
        """
        noise = np.load("noise.npy")
        """
        sig = im[:, 500]
        sig_fft = np.fft.rfft(sig)
        noise_fft = np.zeros_like(sig_fft)
        noise_fft[3] = 50000
        noise_fft[50] = 20000
        noise_fft[200] = 5000
        noise = np.fft.irfft(noise_fft)
        noise /= noise.max()

        nenv = noise_level*signal_level*np.exp(-i/40)
        n = nenv*(2*np.random.random(size=(len(i))) - 1)

        im = (im.T*noise.real).T
        im /= im.max()

        im += np.random.random(size=im.shape)*im.std()/3

        return im

    def _convolve_streak_slit(self, im, t):
        """Blur in the time-domain to account for the width of the streak
        camera

        Parameters
        -----------
        im : 2d np array
            generated sweep
        t : np array
            array corresponding to the time of the sweep
        """
        print("...Convolving streak slit")
        dt = np.diff(t).mean()
        tpx = int(self.tau//dt)

        window = np.ones(shape=tpx)
        return convolve1d(im, window, axis=1)

    def output(self, ray, target, noise=False):
        """Generate the simulated data

        Parameters
        ----------
        ray : Ray class
            the input ray
        target :Target class
            target containing the velocity profile
        noise : bool (optional)
            add in detector noise to the generated image
        """
        I = self._interfear_ray(ray)
        I = self._convolve_streak_slit(I, ray.pulse_length)
        if noise:
            I = self._add_noise(I, ray, target)
        return I


def spatial_var_step(a: "angle", t: "ns", y: "um", max_velocity: "um/ns" = 1):
    """A velocity step-profile which varies linearly in space

    Parameters
    ----------
    a : float
        the slope of the spatially varying profile
    t : float
        the time (in ns) at which to evaluate the velocity
    y : float
        the spatial location at which to evaluate the velocity
    max_velocity : float
        the maximum velocity of the shock

    Returns
    -------
    the velocity determined by the argument parameters
    """
    assert t.shape == y.shape
    v = np.zeros_like(t)
    v[t > -y/a + 3] = max_velocity
    return v


def sin_step(freq, amp, t: "ns", y: "um", max_velocity: "um/ns" = 1):
    """A sinusoidally varying velocity profile in space

    Parameters
    ----------
    freq : float
        the frequency of the spatially varying profile
    amp : float
        the amplitude of oscillations
    t : float
        the time (in ns) at which to evaluate the velocity
    y : float
        the spatial location at which to evaluate the velocity
    max_velocity : float
        the maximum velocity of the shock

    Returns
    -------
    the velocity determined by the argument parameters
    """
    v = np.zeros_like(t)
    v[t > -amp*np.sin(freq*y/(2*np.pi)) + 3] = max_velocity
    return v


def reference_shot(save=False, noise=False):
    """Generate a reference image

    Parameters
    ----
    save : bool (optional)
        save the generated image
    noise : bool (optional)
        add in detector noise

    Returns
    -------
    Pil Image instance
    """
    stationary_target = Target(velocity_equation="stationary")
    ray = Ray(pulse_length=10)
    ray = stationary_target.reflect_off_target(ray)
    etalon = Etalon(1, 1.5195)  # VPF doesn't matter
    interferometer = Interferometer(etalon=etalon)
    ref = interferometer.output(ray, stationary_target, noise)

    ref *= 256/ref.max()
    ref = ref.astype(np.uint8)
    refim = Image.fromarray(ref, mode="L")

    plt.figure()
    plt.imshow(refim, aspect='auto', cmap="gray", extent=(0, 10, -2, 2))
    plt.xlabel("Time [ns]")

    if save:
        refim.save("~/Desktop/ref.jpg", "JPEG")

    return ref


if __name__ == "__main__":
    plt.close("all")

    velocity_equation = lambda t, y: sin_step(20, .5, t, y, max_velocity=1)
    etalon = Etalon(1, 1.5195)
    etalon.set_VPF(2., lambda0=.532)

    # target = Target(velocity_equation="step")
    target = Target(velocity_equation=velocity_equation)
    ray = Ray(pulse_length=10)
    ray = target.reflect_off_target(ray)
    interferometer = Interferometer(etalon=etalon)
    sweep = interferometer.output(ray, target, noise=False)

    plt.figure()
    plt.imshow(sweep, aspect='auto', cmap="gray", extent=(0, 10, -2, 2))
    plt.xlabel("Time [ns]")

    sweep *= 256/sweep.max()
    sweep = sweep.astype(np.uint8)
    im = Image.fromarray(sweep, mode="L")
