import os
import pickle
import json

import matplotlib.pyplot as plt

import visar_core_analysis as core
from etalons import json_to_etalon


def image_from_selector_params(data, params):
    image = core.ImageData(
        data.data[params["str_ij"][1]:params["end_ij"][1],
                  params["str_ij"][0]:params["end_ij"][0]]
    )
    image.set_dimensions(
        xmin=params["str_xy"][0],
        xmax=params["end_xy"][0],
        ymin=params["str_xy"][1],
        ymax=params["end_xy"][1]
    )
    return image


def automated_anaysis(params, raw_data=None):
    if raw_data is None:
        print("...Loading data")
        streak_dims = params["streak_dims"]
        raw_data = core.load_image(
            params["file"],
            xmin=streak_dims["xmin"],
            xmax=streak_dims["xmax"],
            ymin=streak_dims["ymin"],
            ymax=streak_dims["ymax"],
            show=False
        )

    print("...Cropping data")
    cropped_data = image_from_selector_params(
        raw_data, params["data_selection"]
    )

    print("...Generating Spectrogram")
    spec = core.SpectrogramData()
    spec.apply_fft(cropped_data)

    print("...Getting previously selected reference fringes")
    ref_data = image_from_selector_params(
        raw_data, params["ref_fringe_selection"]
    )

    print("...Selecting frequencies to filter")
    ref_spec = core.SpectrogramData()
    ref_spec.apply_fft(ref_data)
    sf, ef = params["freq_selection"]
    filtered_spec, str_freq, end_freq = core.filter_spectrogram(
        spec, str_frequency=sf, end_frequency=ef
    )

    print("...Performing inverse fft")
    filtered_img, filtered_iimg = core.do_inv_fft(filtered_spec, cropped_data)
    print("...Extracting phase")
    phase = core.get_phase(filtered_img, filtered_iimg)
    print("...Generating raw velocity map")
    raw_velocity_map = core.get_velocity_map(phase, ref_data, params["VPF"])

    if params["background_file"]:
        print("...Subtracting background")
        bg_params = params.copy()
        bg_params["file"] = params["background_file"]
        bg_params["background_file"] = None
        bg_velocity_map = automated_anaysis(bg_params)
        raw_velocity_map._data -= bg_velocity_map._data
        velocity_map = raw_velocity_map
    else:
        velocity_map = raw_velocity_map

    print("...Finished")
    return velocity_map


def save_parameters(params, name=None):
    if name is None:
        name = os.path.splitext(params["file"])[0]

    with open('{}.p'.format(name), 'wb') as fp:
        pickle.dump(params, fp)


def load_parameters(path):
    with open(path, 'rb') as fp:
        params = pickle.load(fp)
    return params


def analyze():
    plt.close('all')

    ANALYSIS_PARAMS = {
        "file": None,
        "streak_dims": (),
        "data_selection": (),
        "ref_fringe_selection": (),
        "freq_selection": (),
        "VPF": None,
        "background_file": None,
        "lineout_selection": ()
    }

    with open("streak_camera_conf.json") as f:
        j = json.load(f)
        streak_cam = core.StreakCamera(j["sweep_speed"], j["slit_width"])
        del j

    with open("etalon_conf.json") as f:
        etalons = json.load(f, object_hook=json_to_etalon)

    # Load the raw data
    print("\nLoad png Streak image")
    path = core.prompt_for_path()
    ANALYSIS_PARAMS["file"] = path

    streak_dims = {"xmin": 0, "xmax": 10, "ymin": -290, "ymax": 290}
    raw_data = core.load_image(path, **streak_dims)
    ANALYSIS_PARAMS["streak_dims"] = streak_dims

    # crop out the data the user wants to work with
    print("Select Working data from raw image.")
    data_selector, refrence_selector = core.select_working_data(raw_data)
    cropped_data = data_selector.data()
    cropped_data.show("Cropped Data")
    ANALYSIS_PARAMS["data_selection"] = data_selector.dump_params()

    # Make a spectrogram of the data of the data and apply the filter
    spec = core.SpectrogramData()
    spec.apply_fft(cropped_data)
    spec.show("Spectrogram")

    ref_data = refrence_selector.data()
    ANALYSIS_PARAMS["ref_fringe_selection"] = refrence_selector.dump_params()

    ref_spec = core.SpectrogramData()
    ref_spec.apply_fft(ref_data)

    filtered_spec, str_freq, end_freq = core.filter_spectrogram(spec, ref_spec)
    filtered_spec.show('Filtered Spectrogram')
    ANALYSIS_PARAMS["freq_selection"] = [str_freq, end_freq]

    # Invert the FFT to create the filtered image
    filtered_img, filtered_iimg = core.do_inv_fft(filtered_spec, cropped_data)
    filtered_img.show('Filtered')

    # extract the phase from the arctan of the real and imaginary components
    # to generate the phase image
    phase = core.get_phase(filtered_img, filtered_iimg)
    phase.show("Wrapped Phase")

    # get the etalon used from the user
    vpf = core.vpf_from_etalon(etalons)
    ANALYSIS_PARAMS["VPF"] = vpf
    # vpf = vpf_thru_window(vpf)

    # follow each spacial component to generate a velocity map.
    raw_velocity_map = core.get_velocity_map(phase, ref_data, vpf)

    # Subtract background
    if input("subtract backgound: ").lower()[0] == "y":
        path = core.prompt_for_path()
        ANALYSIS_PARAMS["background_file"] = path
        bg_params = ANALYSIS_PARAMS.copy()
        bg_params["file"] = ANALYSIS_PARAMS["background_file"]
        bg_params["background_file"] = None
        bg_velocity_map = automated_anaysis(bg_params)
        raw_velocity_map._data -= bg_velocity_map._data

    # raw_velocity_map = smooth(raw_velocity_map)
    raw_velocity_plot = raw_velocity_map.show(
        title="Velocity Map", color='jet', colorbar=True
    )

    # Get velocity trace
    time, velocity, lineout = core.get_velocity(
        raw_velocity_map, ax=raw_velocity_plot
    )
    core.plot_velocity(time, velocity, picker=True)
    print("test lineout dump params")
    # ANALYSIS_PARAMS["lineout_selection"] = lineout.dump_params()

    if input("Invert?").lower()[0] == "y":
        velocity *= -1
        plt.cla()
        core.plot_velocity(time, velocity, picker=True)

    while True:

        if input("Add 2 pi shift?").lower()[0] == 'y':
            time, velocity = core.add_fringe_shift(time, velocity, vpf)
        else:
            break
    save_parameters(ANALYSIS_PARAMS, "test_run")

    return time, velocity, ANALYSIS_PARAMS


if __name__ == "__main__":
    plt.ion()
    plt.close('all')
