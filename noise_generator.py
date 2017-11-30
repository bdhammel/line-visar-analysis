import numpy as np
import matplotlib.pyplot as plt
from simulate_data import reference_shot
from scipy.interpolate import interp1d


def add_freq(ref_fft, freq, amp):
    norm = ref_fft.sum()
    ref_fft[freq] = amp
    ref_fft *= (ref_fft.sum()/norm)
    return ref_fft

def generate_noise(ref_fft):
    noise = np.fft.irfft(ref_fft)
    return noise

def mult_sig(sigA, sigA_w, sigA_b, sigB):
    norm = sigB.sum()
    sig = (sigA_w * sigA.real + sigA_b) * sigB
    return sig * norm/sig.sum()

plt.close('all')

ref_im = reference_shot(save=False, noise=False)
ref_sig = ref_im[:,500]
plt.figure("fundamental sig")
plt.plot(ref_sig)
plt.title("fundamental sig")
plt.figure()
fundamental_fft = np.fft.rfft(ref_sig)
plt.plot(fundamental_fft.real)
ref_fft = np.copy(fundamental_fft)
noise_fft = np.zeros_like(fundamental_fft)

real_im = plt.imread("/Users/bdhammel/Desktop/2017060202_ref.png")
_real_sig = real_im[:,400:500].mean(axis=1)
_x = np.arange(len(_real_sig))
_real_sig_f = interp1d(_x, _real_sig)
real_sig = _real_sig_f(np.linspace(0, _x[-1], len(ref_sig)))
plt.figure()
plt.plot(real_sig)
real_fft = np.fft.rfft(real_sig)


noise_level = .8
sig_level = ref_sig[5:].max()

nenv = noise_level*sig_level*np.exp(-i/40)
n = nenv*(2*np.random.random(size=(len(i)))- 1)

im = (im.T * noise.real).T
im/=im.max()

im += np.random.random(size=im.shape)*im.std()/3

    
noise = np.fft.irfft(noise_fft)
s = mult_sig(noise, .75, 50, ref_sig)
plt.plot(s/s.max()*.45)
