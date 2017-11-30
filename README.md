# VISAR Scripts

For a discussion of the theory, please visit [my blog post](https://bdhammel.github.io/2017/06/10/line-visar-analysis.html)

## Analysis

Run the analysis script in an iPython terminal. Click-and-drag the image (png format) in the terminal.

~~~python
>>> %run visar_analysis.py
>>> time, velocity, analysis_params = analyze()

Load png Streak image
Drag and drop location into command prompt
path > /<Path>/smith.png
~~~

### Select working data and reference

This will open up the PNG image in a matplotlib figure. The figure has event handeling enabled, and you will be propted to `Select Working data from raw image` by drawing a rectangle (red) on the figure by clicking and dragging. You must drag from the bottom-left to top-right.

![](/Users/bdhammel/Documents/research_programing/visar_scripts/media/raw_data.png)

Press <kb>enter</kb> in the terminal to select the region of data. You will then be propted to select a region of reference fringes, drag away from the start of the rectangle (green). Press <kb>enter</kb> in the terminal to confirm selection.

### Select the fundamental frequency

A new window will open for you to select the fringe-spacing frequency of interest. Click on a location of the plot to select the minimum frequency (green), press enter (at the window focus, not the terminal), and select the maximum frequency (red). Press enter in the terminal to end selections.

![](/Users/bdhammel/Documents/research_programing/visar_scripts/media/freq.png)

### Pick Etalon

A menu will prompt you to choose the VPF. The VPFs listed in this menu are determind my the JSON config file `etalon_conf.json`.

~~~bash
Pick an etalon thickness:

	[0]____________1e+03 mm________46.29 km s-1
	[1]____________2e+03 mm________23.14 km s-1
	[2]____________5e+03 mm________9.258 km s-1
	[3]________1.144e+04 mm________4.047 km s-1
	[4]__________1.5e+04 mm________3.086 km s-1
	[5]__________2.5e+04 mm________1.852 km s-1
	[6]__________4.5e+04 mm________1.029 km s-1
	[7]_________9.97e+03 mm________5.175 km s-1
	[8]_________2.53e+04 mm________2.039 km s-1
	[9]______________custom
> 9
Enter custom vpf: 5.4603
~~~

### Subtract Background (optional)

Enter "yes" in the terminal to subtract background fringes. You will be prompted to drag-and-drop a reference image. This will then be processed with the same selections you made to the data. (These parameters are stored in the `ANALYSIS_PARAMS` dictionary. 

### Take a line-out of the data

Using the click-and-drag technique, select a rectangle of interest (green) from the veloctiy map.

![](/Users/bdhammel/Documents/research_programing/visar_scripts/media/vmap.png)

`analyse()` will now return time and velocity arrays from your selection, along with a dictionary of the parameters you selected during the analysis

### EXTRA

You will then be propted to add fringe shifts, or invert the velocity trace


## Recalling an old Analysis


### Saving analysis parameters

At the end of the analysis, a `analysis_param` dictionary will be returned. Save this dictionary using the method call 

~~~python
save_parameters(analysis_params, name="name/of/file")
~~~

This will generat a pickle file which can be recalled later

### Loading alaysis parameters

To recall an old analysis session, run the automation method with the pickle analysis_params passed as an argument 

~~~python
automated_analysis(params="name/of/file.npy")
~~~
