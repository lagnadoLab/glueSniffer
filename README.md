# glueSniffer
Glue Sniffer Analysis

# Load Movies
The load movies button creates a popup listbox allowing the user to select a linescan movie to analyze. The current version of the code only allows for linear linescans (ie. no curves in the linescan). 

# Enter Params
This button creates a popup box allowing the user to enter the parameters of the Wiener Filter - tau rise and tau fall. Importantly, these parameters must be initialized with the button before deconvolution is completed. 

# Split Channels
This button allows the user to split the channels of the linescan video. In our recording setup, stimulus frames are interleaved with recording frames. Splitting the channels seperates these two aspects.

# Linescan
This button converts the linescan video into a pixels x time matrix for easier analysis. 

# New Define AZs
Creates GUI for defining active zones. A temporal average of the linescan matrix is displayed, and users place cursors on areas where they believe ROIs are centered. Igor's curve-fitting software is then used to fit the temporal average to a sum of Gaussians, where means are defined by cursor placement (with some margin for error). This can be used multiple times if a decomposition doesn't work.

# Old Define AZs
Old version of defining AZs. As in the new version, a temporal average is displayed. In this version however, users simply place cursrs to define the ROIs. The ROI is then taken to be the area between each pair of cursors.

# New Temporal Profile
Creates a temporal profile for each ROI in New Define ROIs by weighted averaging based on the Gaussian decomposition. Thus, values nearer the center of each ROI are weighted heavier than outside areas. Calculates a 'deconvolved' wave based on each ROIs temporal profile and the Wiener filter with parameters initialized from Enter Params.

# Old Temporal Profile
Creates a temporal profile by uniformly averaging all points within the ROI defined in Old Define AZs. Deconvolves in same fashion as New Temporal Profile

# Threshold
Initializes an event detection GUI. Events are defined as local maxima of the deconvolved wave above a certain user-set threshold. The GUI graphically displays detected events based on a threshold set with a slider bar.
