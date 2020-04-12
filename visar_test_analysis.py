"""
Things to vary:
    target:
    - fringe jumps
    - spatial variation
    - number of fringes

    analysis:
    - frequency range to select
    - index of refraction of window
    - slit width

    Noise
"""
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

from visar_analysis import automated_anaysis, load_parameters
from simulate_data import (Ray, Etalon, Target, Interferometer,
                           spatial_var_step)


if __name__ == "__main__":

    velocity_equation = lambda t, y: spatial_var_step(3, t, y, max_velocity=3)
    # velocity_equation = lambda t, y : sin_step(20, .5, t, y, max_velocity=1)

    plt.close('all')

    EPOCHS = range(1)

    avg = []
    std = []

    print("...Initializing")
    target = Target(velocity_equation="step")
    ray = Ray(pulse_length=10)
    etalon = Etalon(1, 1.5195)
    etalon.set_VPF(2., lambda0=.532)
    interferometer = Interferometer(etalon=etalon)

    print("...Calculating phase shift")
    ray = target.reflect_off_target(ray)

    sweep = interferometer.output(ray, target)

    sweep *= 256/sweep.max()
    sweep = sweep.astype(np.uint8)
    im = Image.fromarray(sweep, mode="L")
    im.save("/Users/bdhammel/Desktop/test.jpg", "JPEG")

    p = load_parameters("/Users/bdhammel/Desktop/test.p")
    p['background_file'] = "/Users/bdhammel/Desktop/ref.jpg"
    vmap = automated_anaysis(p)

    tm, _ = vmap.physical_to_pixel(6, 0)
    vmax = vmap.data[800:950, tm:]

    avg.append(vmax.mean())
    std.append(vmax.std())
