import numpy as np
from scipy.signal import savgol_filter  
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os

"""Line VISAR analysis script 

__author__ = "Ben Hammel"
__email__ = "bdhammel@gmail.com"
__status__ = "Development"

References
-----------
[1] P. M. Celliers, D. K. Bradley, G. W. Collins, D. G. Hicks, T. R. Boehly, 
and W. J. Armstrong, "Line- imaging velocimeter for shock diagnostics at the 
OMEGA laser facility," Review of Scientific Instruments, vol. 75, pp. 4916–4929, 
nov 2004
"""

class StreakCamera:
    def __init__(self, sweep_window, slit_width, slit_opening=None):
        self.sweep_window = sweep_window
        self.slit_width = slit_width


class ImageData:
    """Holder class for spatial image data

    Attributes
    ----------
    data (2D numpy array): image data
    xx (2D nparray) : x values for each pixel in image data
    yy : y values for each pixel in image data
    """

    def __init__(self, img):
        """Load image data and construct x and y values for each pixel 

        Args
        ----
        img (2D list) : 2d array of image data

        """
        self._data = np.asarray(img)


    def set_dimensions(self, xmin, xmax, ymin, ymax):
        """set the physical dimensions of the image
        Args
        ----
        xmax (int): the length of the sweep window 
        ymax (int): the size of the slit on the streak camera as it relates to 
            the image focal size
        """
        jx, ix = self.data.shape
        self.x = np.linspace(xmin, xmax, ix)
        self.y = np.linspace(ymin, ymax, jx)

    def copy_dimensions(self, data):
        """Copy dimensions from another ImageData object
        """
        self.x = data.x
        self.y = data.y

    @property
    def data(self):
        """raw image data
        """
        return self._data

    def pixel_to_physical(self, px, py):
        """Convert a pixel location to physical location
        """
        xx, yy = np.meshgrid(self.x, self.y)
        return xx[px, py], yy[px, py] 

    def physical_to_pixel(self, x, y):
        """Convert physical location to a pixel location
        """
        px = np.argmin(np.abs(self.x-x))
        py = np.argmin(np.abs(self.y-y))
        return px, py

    @property
    def aspect(self):
        """Aspect ratio of the image
        Ration of x pixels / y pixels 

        Returns
        -------
        float : xpixels/ypixels
        """
        return self._data.shape[0]/self._data.shape[1]

                
    def show(self, title=None, color='gray', frame=True, figsize=(6,4), colorbar=False):
        """Display the image

        Args
        ----
        title (str) : Title of the figure window
        color (str) : Color map to use
        frame (bool) : display ticks and axis labels
        figsize (int, int) : aspect ratio of the figure 
        """
        
        fig = plt.figure(title, figsize=figsize)
        fig.set_tight_layout(True)
        ax = fig.add_subplot(111)

        if frame:
            ax.set_xlabel("Time [ns]")
            ax.set_ylabel("X [um]")
            extent = (self.x.min(), self.x.max(), self.y.min(), self.y.max())
        else:
            ax.set_axis_off()
            extent = None

        implot = ax.imshow(self._data, aspect='auto', extent=extent, origin="lower")
        implot.set_cmap(color)

        if colorbar:
            fig.colorbar(implot)

        plt.show()

        return implot.axes

class SpectrogramData(ImageData):
    """Extension of ImageData, specific to displaying a spectrogram 
    """

    def __init__(self, spec=None, wavenumbers=None, x=None, *args, **kwargs):
        if spec is not None:
            self.k = wavenumbers
            self.x = x
            super().__init__(spec, *args, **kwargs)

    def apply_fft(self, img):
        """Create a spectrogram from given data

        Create a spectrogram from image data, save the original x (time) values
        from the original image data object, and create a wave number array 
        from the spatial values (y)

        Args
        ----
        img (ImageData): object containing the data to Fourier transform 
        """
        self._data = np.fft.fft(img.data, axis=0)
        self.k = get_wavenumber(img.y, self._data.shape[0])
        self.x = img.x

    @property
    def pdata(self):
        """Pretty data - rearrange the y values so the frequencies are plotted
        in an intuitive way
        # -N/2 ... N/2 range 
        """
        inds = np.argsort(self.k)
        return np.abs(self.data)[inds]

    @property
    def pwavenum(self):
        """Pretty wavenumber - corresponds to pretty data
        # -N/2 ... N/2 range 
        """
        return sorted(self.k)

    def show(self, title, color='hot', figsize=(6,4), frame=True):

        fig = plt.figure(title, figsize=figsize)
        fig.set_tight_layout(True)
        ax = fig.add_subplot(111)

        if frame:
            ax.set_xlabel("Time [ns]")
            ax.set_ylabel("K [mm${}^{-1}$]")
            extent = (self.x.min(), self.x.max(), self.k.min(), self.k.max())
        else:
            ax.set_axis_off()
            extent = None

        implot = ax.imshow(self.pdata, aspect='auto', extent=extent, origin="lower")
        implot.set_cmap(color)
        plt.ylim(-0.1, 0.1)

        plt.draw()

        return implot.axes


class IndexSelector:
    """
    Return the index of the x value clicked on the array

    Attributes
    ----------
    _ax (Axis): axis the plot is drawn on
    xdata (nparray): xvalue of the line on the graph
    _cursor_index (int): the index value the cursor was placed at
    _line : object drawn on plot after user click
    _color: color the line should be

    Example
    -------
    line, = ax.plot(xs, ys, 'o', picker=5) # setting picker is important!
    sel = IndexSelector(ax, line=line)

    """
    def __init__(self, line, color='g'):
        """
        Args
        ----
        line (matplotlib line object) :
        color (str) :
        """
        self._ax = line.axes
        self.xdata = line.get_xdata()

        self._cursor_index = None
        self._line = None
        self._color = color

        self.cidclick = self._ax.figure.canvas.mpl_connect("pick_event", self._on_pick)
        self.cidkey = self._ax.figure.canvas.mpl_connect("key_press_event", self._key_press)

    def _on_pick(self, event):
        """When there is a mouse click, record the location of the click, and draw
            a line. 
        """
        try:
            self._cursor_index = event.ind[0]
        except ValueError:
            pass
        else:
            self.draw(self.xdata[self._cursor_index], color=self._color)

    def _key_press(self, event):
        sys.stdout.flush()
        if event.key == 'return' or 'enter':
            self.disconnect()

    def remove_line(self):
        """Remove the vertical line if it's there
        """
        if self._line:
            self._line.remove()

    def draw(self, location, color='g'):
        """Draw a line on the graph
        """
        self.remove_line()
        self._line = self._ax.axvline(location, color=color)
        plt.draw()

    def disconnect(self):
        """Return true of successful disconnect
        """
        self._ax.figure.canvas.mpl_disconnect(self.cidclick)
        self._ax.figure.canvas.mpl_disconnect(self.cidkey)

    def get_index(self):
        """return the index of the line
        """
        return self._cursor_index

    def get_xdata(self):
        """return the xdata of the pick event 
        """
        return self.xdata[self.get_index()]


class DataSelector1D(IndexSelector):
    """Select data between two cursor positions
    returns the starting and ending index on a set of that a
    """

    def __init__(self, line):
        self.line = line
        self.starting_cursor = True
        super().__init__(self.line, color='g')

    def _key_press(self, event):
        """if enter is pressed, move to next cursor.
        Disconnect event handling if ending cursor is selected
        """
        if self.starting_cursor is True:
            self.starting_index = self.get_index()
            self.disconnect()
            super().__init__(self.line, color='r')
            self.starting_cursor = False
        else:
            self.ending_index = self.get_index()
            self.disconnect()

    def get_indices(self):
        """return the starting and ending indices to the user
        """
        return self.starting_index, self.ending_index

    def get_xdata_pts(self):
        """Return the x data points at the selected locations
        Returns
        -------
        (float, float) : x1, x2 data points at locations of cursor 1 and cursor 2
        """
        idx = self.get_indices()
        return self.xdata[list(idx)]


class DataSelector2D:
    """Select a 2d section of data on an image

    Allow user to draw a rectangle on the image
    information enclosed in the rectangle is returned to the user
    
    Rectangle is drawn by left clicking the mouse, and moving,
    releasing the left click sets the rectangle.
    when enter is pressed, disconnect the canvas

    Attributes
    ----------
    _ax (pyplot axis): Axis the image is being displayed on
    _draw (boolean): Should the selector track mouse movement to draw the 
        rectangle in real time
    _rec (pyplot Rectangle): Artist object drawn on figure
    _color (str): color of rectangle g, r, b, etc...
    str_xy (tuple): starting coordinates of rectangle (upper left)
    end_xy (tuple): ending coordinates of rectangle (bottom right)
    width (int): width of rectangle
    height (int): height of rectangle
    """

    def __init__(self, img, ax=None, color='r'):
        if not ax:
            self._ax = plt.gca()
        else:
            self._ax = ax

        self._draw = False
        self._rec = plt.Rectangle((0,0),0,0)
        self._color = color
        self.img = img

        # add initial rectangle to canvas
        self._ax.add_patch(self._rec)
        
        # initialize key-bindings and event handling
        self.cidpress = self._ax.figure.canvas.mpl_connect("button_press_event", self._mouse_press)
        self.cidrelease = self._ax.figure.canvas.mpl_connect("button_release_event", self._mouse_release)
        self.cidmove = self._ax.figure.canvas.mpl_connect("motion_notify_event", self._mouse_move)

    def is_valid(self):
        """Check if a valid selection has been made
        Returns
        -------
        (bool) : True if the height and width are positive, False otherwise
        """
        return (self.height > 0) and (self.width > 0)

    def disconnect(self):
        """Disconnect the canvas when "enter" is pressed
         set "waiting_for_entry" to False:
           used in main loop, so that the code and move to to next step
        def _key_press(self, event, force_disconnect=False):
        """
        if self.is_valid():
            self._ax.figure.canvas.mpl_disconnect(self.cidpress)
            self._ax.figure.canvas.mpl_disconnect(self.cidrelease)
            self._ax.figure.canvas.mpl_disconnect(self.cidmove)
            return True
        else:
            print("Selection must start in lower left and go to upper right")
            return False

    def _mouse_press(self, event):
        """Activate on mouse press event

        Record the starting place of the rectangle, enable the rectangle to be 'drawn'
        """
        self.str_xy = np.array([event.xdata, event.ydata])
        self._draw = True

    def _mouse_move(self, event):
        """Event as the mouse is moved across the image
        active only after initial mouse press event

        Draw the rectangle as the mouse is moved and record the location
        """
        if self._draw and event.inaxes:
            self.end_xy = np.array([event.xdata, event.ydata])
            self._draw_rec()

    def _draw_rec(self):
        """Render the Rectangle to the figure
        """
        # Remove existing rectangle
        if self._rec is not None:
            self._rec.remove()

        # generate new rectangle corresponding to new location of the mouse
        self._rec = plt.Rectangle(
            self.str_xy, self.width, self.height, 
            fill=False, linewidth=1.5, color=self._color)

        # if the user is not selecting from the bottom left to upper right, 
        # give visual to show invalid selection 
        if not self.is_valid():
            self._rec.fill = True

        self._ax.add_patch(self._rec)
        plt.draw()

    def _mouse_release(self, event):
        """Run on a mouse release event
        disable drawing, the user has chosen the rectangle
        """
        self._draw = False

    @property
    def height(self):
        """The height of the current selected region

        Calculate the height based on the starting location and ending location
        of the mouse

        Returns
        -------
        float
        """
        return (self.end_xy - self.str_xy)[1]

    @property
    def width(self):
        """The width of the current selected region

        Calculate the width based on the starting location and ending location
        of the mouse

        Returns
        -------
        float
        """
        return (self.end_xy - self.str_xy)[0]

    def data(self):
        """return the data enclosed by the rectangle

        TODO: Check how xmin and xmax are being passed to the next image data object

        Returns
        -------
        ImageData object : cropped image
        """
        xmin, ymin = self.str_xy
        xmax, ymax = self.end_xy

        # Get the indexes of the starting and ending positions
        self.str_ij = self.img.physical_to_pixel(*self.str_xy)
        self.end_ij = self.img.physical_to_pixel(*self.end_xy)

        # get the image back from the axis
        working_image = ImageData(
                self.img.data[
                    self.str_ij[1]:self.end_ij[1], 
                    self.str_ij[0]:self.end_ij[0]
                    ])
        working_image.set_dimensions(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

        return working_image

    def get_ax(self):
        """Return the current axes being operated on
        """
        return self._ax

    def dump_params(self):
        """Dump out parameters of selector
        """
        return {"str_ij":self.str_ij, 
                "end_ij":self.end_ij,
                "str_xy":self.str_xy, 
                "end_xy":self.end_xy}

class LockedSelector(DataSelector2D):
    """Select a section of data with specific constraints

    Allows user to draw a rectangle on the data image under specific constraints
    i.e. only draw a rectangle starting at (x,y)
         only draw a rectangle of width (or height) ____
    """

    def __init__(self, img, lock_start=None, lock_width=None, lock_height=None,
            **kwargs):
        super().__init__(img, color='g', **kwargs)

        self.lock_start = lock_start
        self.lock_width = lock_width
        self.lock_height = lock_height
                        
    def _mouse_press(self, event):
        """Initiate the drawing of the selector

        If a start location was passed to the instance during initialisation, 
        use that. Otherwise, use the mouse click location.
        """
        super()._mouse_press(event)

        if self.lock_start is not None:
            self.str_xy = self.lock_start

    def _mouse_release(self, event):
        """Run on a mouse release event
        disable drawing, the user has chosen the rectangle

        Update end location to be consistent with the locked dimensions
        """
        super()._mouse_release(event)
        self.end_xy = self.str_xy + np.array([self.width, self.height])

    @property
    def width(self):
        """The width of the current selected region

        Returns
        -------
        float : if width is specified during initialization, return this with, 
            otherwise call the DataSelector2D width method 
        """
        if self.lock_width is not None:
           return self.lock_width
        else:
            return super().width

    @property
    def height(self):
        """The height of the current selected region

        Returns
        -------
        float : if height is specified during initialization, return this with, 
            otherwise call the DataSelector2D height method 
        """
        if self.lock_height is not None:
            return self.lock_height
        else:
            return super().height

    def force_return(self):
        """return a selection without user input 
        """
        self.str_xy = self.lock_start
        self.end_xy = self.str_xy + np.array([self.lock_width, self.lock_height])

        super().disconnect()
        return super().data()


def get_wavenumber(y, N):
    width = y.max() - y.min()
    width = np.linspace(0, width, N)
    step_size = width[1]-width[0]

    # hz = [0, 1, ..., N/2-1, -N/2, ..., -1] / (d*N)         if N is even
    # hz = [0, 1, ..., (N-1)/2, -(N-1)/2, ..., -1] / (d*N)   if N is odd
    return np.fft.fftfreq(N, d=step_size)


def prompt_for_path(path=None):
    """Ask user to submit path to the desired directory or file

    Prompts user to drag and drop file or directory into the terminal window, this
        path is 'cleaned' and returned so that the desired object can be found

    Args
    ----
    path (str, optional): The path can be submitted via the function call, 
            and the user wont be prompted through the terminal 

    Returns
    -------
    (str): path the user submitted, with the correct use of separators, and any
    leading, or trailing, spaces removed. 
    """

    if not path:
        print("Drag and drop location into command prompt")
        path = input("path > ")

    print('')  # new line for spacing

    if os.name is 'posix':
        return path.strip().replace('\\', '')
    else:
        print("Might not work on Windows")
        return path 


def wait_for_input(selector):
    """Hold the program until the user has made a valid selection

    Check if data selector returns a valid disconnect, if not, continue to wait

    Args
    ----
    (DataSelector) : Data selector object
    """
    while True:
        input("...waiting ")
        if selector.disconnect():
            break

def load_image(path, xmin, xmax, ymin, ymax, show=True):
    """Load the Raw Data

    Load a png image into a numpy array
    image must be uploaded such that fringes are running up and down

    Parameters
    ----------
    path (str) : cleaned path loaded from user

    Return
    ------
    ImageData Object

    """
    img_array = plt.imread(path)

    img = ImageData(img_array)
    img.set_dimensions(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    if show:
        img.show('Raw Data', figsize=(6, 6/img.aspect))

    return img


def select_working_data(img):
    """Return the section of data the user wants to work with
 
    Return full section of working data and sub section for constant fringes
    to be used as a reference in FFT analysis
    """
    data_selector = DataSelector2D(img)

    # hold program until the user has selected the data and pressed "enter"
    wait_for_input(data_selector)

    print("Select constant fringes for FFT reference")
    refrence_selector = LockedSelector(
            img, 
            lock_start=data_selector.str_xy,
            lock_height=data_selector.height)
    wait_for_input(refrence_selector)
    
    return data_selector, refrence_selector

def show_averaged_fft(ref_spec):
    """Display a section of a spectrogram, averaged across each spectral dimension 

    Args
    ----
    ref_spec (SpectrogramData) : Reference spectrogram of constant fringes 

    Returns
    -------
    (matplotlib line) : line object of the averaged fft, with a "picker" value set
    """
    plt.figure("Reference frequencies", figsize=(5,3))
    avg_ref_spec = np.average(ref_spec.pdata, axis=1)
    l, = plt.plot(ref_spec.pwavenum, avg_ref_spec, picker=5)
    plt.xlabel("k")
    plt.xlim(0, 0.2)
    plt.draw()
    plt.tight_layout()
    return l

def filter_spectrogram(spec, ref_spec=None, str_frequency=None, end_frequency=None):
    """Remove unwanted frequencies by setting values in the spectrogram to 0

    Args
    ----
    spec (SpectrogramData) :  Spectrogram of the selected working data
    ref_spec (SpectrogramData) : Reference spectrogram of constant fringes 
    start_frequency (float): the starting frequency to keep (passed during analysis of background)
    end_frequency (float): the last frequency to keep (passed during analysis of background)
    """

    # Takes the first part of the image (constant fringes) and averages the FFT
    # lets the user select the starting and ending frequencies they want to keep
    if ref_spec is not None:
        fft_avg_line = show_averaged_fft(ref_spec)

        print("Select frequencies of interest")
        index_selector = DataSelector1D(fft_avg_line)
        input("...waiting ")
        index_selector.disconnect()

        str_frequency, end_frequency = index_selector.get_xdata_pts() 


    # Get all of the indices of the wave number array between the points selected
    # by the user
    filtered_spec = np.copy(spec.data)
    filter_indices = np.where((spec.k < str_frequency) | (spec.k > end_frequency))
    filtered_spec[filter_indices,:] = 0

    filtered_spec = SpectrogramData(
            filtered_spec, 
            wavenumbers=spec.k,
            x=spec.x)

    return filtered_spec, str_frequency, end_frequency

def do_inv_fft(spec, cropped_data):
    """Preform an inverse fft of each row of the spectrogram
    recreate the image with the filter applied

    Args
    ----
    spec (SpectrogramData) : Filtered spectrogram
    cropped_data (ImageData) : selection of the data chosen by the user

    Returns
    -------
    filtered_img (ImageData) : filter image from real numbers
    filtered_iimg (ImageData) : filter image from imaginary numbers
    """

    # both real and imaginary parts are needed to do the phase unwrap
    filtered_iimg = np.fft.ifft(spec.data, axis=0).imag
    filtered_img  = np.fft.ifft(spec.data, axis=0).real

    filtered_img  = ImageData(filtered_img) 
    filtered_img.copy_dimensions(cropped_data)

    filtered_iimg = ImageData(filtered_iimg)
    filtered_iimg.copy_dimensions(cropped_data)

    return filtered_img, filtered_iimg

def get_phase(img, iimg):
    """Get the wrapped phase 
 
    W(x) = archtan( sin(x) / cos(x) )

    Args
    ----
    img  (ImageData) : is the real image
    iimg  (ImageData) : is the imaginary image 

    Returns
    -------
    (ImageData) : wrapped phase image
    """
    phase = np.arctan(iimg.data/img.data)
    phase = ImageData(phase)
    phase.copy_dimensions(img)

    return phase 

def vpf_from_etalon(etalons):
    """Calculate the velocity per fringe 
    eqn (2) from [1]
    """
    print("Pick an etalon thickness:\n")
    for i, etalon in enumerate(etalons):
        print("\t[{}]{thickness:_>20}{vpf:_>20}".format(
            i,
            thickness="{:.4} mm".format(etalon.thickness),
            vpf="{:.4} km s-1".format(etalon.vpf),
            ))
    print("\t[{}]{:_>20}".format(i+1, "custom"))

    choice = int(input("> "))
    try:
        vpf = etalons[choice].vpf
    except:
        vpf = float(input("Enter custom vpf: "))

    return vpf

def vpf_thru_unshocked_window(vpf):
    """Correct the VPF to a window'd target

    Only applicable to measuring a reflective shock front

    Equation 3

    """
    windows = [
            {"material":"None", "correction":0},
            {"material":"LiF", "correction":0},
            {"material":"CH", "correction":1.57},
            ]

    print("Window material? \n")
    for i, window in enumerate(windows):
        print("\t[{}]{material:_>20}".format(
            i, material=window["material"]))
    print("\t[{}]{:_>20}".format(i+1, "custom"))

    choice = int(input("> "))
    vpf = vpf/(1+windows[choice]["correction"])

    return vpf

def get_velocity_map(img, ref, vpf):
    """Construct velocity map from phase map using VPF

    Velocity is multiplied by -1 to be consistent with plotted images - images
    are flipped, as the origin is set to be at the bottom left during plotting

    TODO:
    Do not set diff to 0, set it to interpolate between the previous and future
    values...

    Args
    ----
    img (ImageData) : 
    ref ()
    vpf (float) : velocity per fringe shift - dependent on the etalon chosen 

    Returns
    -------
    """
    avg_to = len(ref.x)

    wrapped_phase = np.copy(img.data)

    # find the phase change of each spatial point relative to the reference 
    # phase, at that same spatial point
    rel_phase = wrapped_phase \
                - np.mean(wrapped_phase[...,:avg_to], axis=1, keepdims=True)

    dif = np.diff(rel_phase, axis=1)
    dif[np.where(dif > 1.5)] -= np.pi
    dif[np.where(dif < -1.5)] += np.pi

    unwrapped_phase = np.cumsum(dif, axis=1)

    # Generate velocity image based on VPF
    velocity_img = ImageData(unwrapped_phase*vpf/(2*np.pi))
    velocity_img.copy_dimensions(img)

    return velocity_img

def smooth(velocity_map):
    """Smooth the image in the spatial direction
    """
    a = savgol_filter(velocity_map.data, window_length=201, polyorder=2, axis=0)
    velocity_map._data = a
    return velocity_map

def get_velocity(velocity_map, ax, lineout_height=5.):
    """Take a lineout of the velocity map 

    let the user select a section of the velocity map.
    average the values for each temporal point, and return a 1d array.:
    construct the time values based on the ratio used in the image properties.

    Args
    ----
    velocity_map (ImageData): 2D phase map of velocity
    lineout_height (negative int): height of lineout used
        to select average of velocity.     
    """
    lineout_selector = LockedSelector(velocity_map, ax=ax, lock_height=lineout_height)

    print("\nSelect a line out. width: {}".format(lineout_height))
    wait_for_input(lineout_selector)
    
    lineout = lineout_selector.data()
    velocity = np.average(lineout.data, axis=0)

    return lineout.x, velocity, lineout 

def plot_velocity(time, velocity, picker=False):
    """Helper function to plot velocity trace

    Args
    ----
    time (nparray float): time values
    velocity (nparray float): velocity values
    picker (bool) : enable picker so that a selector can be used on the plot
    """

    fig =  plt.figure("Velocity")
    fig.clear()

    if picker:
        plt.plot(time, velocity, picker=2)
    else: 
        plt.plot(time, velocity)

    plt.grid()

    plt.xlabel("Time (ns)")
    plt.ylabel(r"Velocity $({\rm km \, s^{-1}})$")
    plt.tight_layout()
    plt.draw()

def add_fringe_shift(time, velocity, vpf):
    """add a full fringe shift to user selected point
    """
    ax = plt.gca()
    line = ax.lines[0]

    # let user select discontinuity
    picker = DataSelector1D(line)
    input("Select both cursors...")

    i1, i2 = picker.get_indices()
    delta_v = velocity[i1]-velocity[i2]

    # add a full fringe shift, and remove some data points before the shift (these
    # would have been generated incorrectly during the analysis)
    print("up or down?")
    choice = input().lower()
    print("Number of jumps?")
    mult_factor = float(input())

    if choice[0] == 'u':
        velocity[i2:] += (mult_factor*vpf) # + delta_v)
    else:
        velocity[i2:] -= (mult_factor*vpf) # + delta_v)

    delete_i = np.arange(i1,i2)
    velocity = np.delete(velocity, delete_i)
    time = np.delete(time, delete_i)

    # remove old velocity trace and replot
    plt.cla()
    plot_velocity(time, velocity, picker = True)

    return time, velocity

def subtract_background():
    print("\nLoad background streak image")
    path = prompt_for_path()
    ANALYSIS_PARAMS["background_file"] = path
    raw_background = load_image(path, *ANALYSIS_PARAMS["streak_dims"])
    cropped_data = None
    spec = SpectrogramData()
    spec.apply_fft(cropped_data)
    ref_data = None
    ref_spec = SpectrogramData()
    ref_spec.apply_fft(ref_data)

