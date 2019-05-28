# moveflow
Compute optical flow + free form deformations to capture the movement of the organism. Align image slices in the image stack. Tracking neuron and quantify the Calcium activity. Edit (divide) traces into neurons. Compute the stress over the image volume based on the given optical flow.

Tracking neuron in calcium image datasets uses Markov Network where dendrites are modeled as the articulated body. MRF is optimized using fusion moves (QPBO + alpha expansion).

Example
================
* trackGFPRFPXYTdata.m example for 3D data (2D spatial XY + time T)
* trackGFPRFPXYZTdata.m example for 4D data (3D spatial XYZ + time T)

Folder Hierarchy
================
* Full_Flow_Source_Code/myfullflow.m. Parameters for each dataset are bundled in switch statement. Produce the optical field in 2D of collapsed image stack by maximum intensity over time in variables `U` and `V`, while collapsed image stack is stored in `Is`.
Compute Optical Flow using [Full FLow](http://web.stanford.edu/~cqf/fullflow/) from 'Full Flow: Optical Flow Estimation By Global Optimization over Regular Grids' by Qifeng Chen and Vladlen Koltun. 

* EditTrace/GUI_edit_trace.m. The GUI for removing edges in the given neuron trace. It comes with the installer, *Edit Trace.mlappinstall*.

* myffd.m. Compute the stress over the image stack based on the first derivative of the optical flow with respect to pixel in X direction (compute partial derivative of the vector field, `U` and `V`, in X direction). Results are stored in videos.

* demon_registration_version_8f/registerstacks.m. Align image slices in the image stack using [demon image registration](https://www.mathworks.com/matlabcentral/fileexchange/21451-multimodality-non-rigid-demon-algorithm-image-registration). Similarly, parameters for each dataset are bundled in switch statement. The code is very slow. Computation is reduced by register only image stack that has movement, which is determined by the total magnitude of the optical flow exceeding the threshold. Then, +/- one frame in case the movement was not detected. The result stored as MATLAB file (.mat) with one cell array variable `Vregis` storing the aligned image stacks only at specific time step.

* pick_frame.m. Display the collapsed image stack over time to help user pick frames to trace neurons for tracking. Selected image stacks are stored in `<neuron_name>_frames` folder.

* track_neuron.m. Tracking neuron by moving the trace along the optical flow. Traces are computed manually at one specific frame with the highest clarity and are stored in the folder `<neuron_name>_frames`. Then, traces are moved along optical flow both forward and backward until they move out of bound. Then, pixels around the trace are used as the filter mask for quantifying the Calcium activity. The result is stored as video.

* mrftrack.m. Main tracking function. It models the trace as the articulated body and solves the transformation using MRF with fusion moves (QPBO + alpha expansion)
```matlab
% MRFTRACK
% By Sarun Gulyanon
%
% Usage
%
% [final_maskA, final_maskB] = main_coseg( IA, IB, WA, WB, maskA, maskB,...
%                               somaA, somaB, maskA_gold, maskB_gold, opt )
%
% This function tracking neuron in calcium image datasets using Markov 
% Network where dendrites are modeled as the articulated body. MRF is 
% optimized using fusion moves (QPBO + alpha expansion).
%
% Inputs:
%   swc           trace at time t_{n-1} in SWC format
%   Iori          calcium image stack at time t_n
%   soma          soma location at time t_n
%   soma_prev     soma location at time t_{n-1}
%   swc_model     trace at time t_0 in SWC format
%   soma_model    soma location at time t_0
%   opt      miscellaneous parameters.
%       'alpha_d'   Tune dynamic model. Default 1/5.
%       'alpha_t'   Tune transition similarity. Default 1/10.
%       'DEBUG'     Enable debugging mode. Default false.
%       'MAXDIST'   Maximum displacement. Default 10.
%
% Outputs:
%   swc           Trace at time t_n in SWC format
%   flag          True if the neuron moves out-of-bound, otherwise False.
```

Video Results
================
Videos show the image-processing results of 4D neuron data that captures the locomotion of the organism. Each data sequence has 6 files (if applicable). Names of videos follow this pattern: 
* `<neuron_name>_fullflow.avi` shows the optical flow at every frame.
* `<neuron_name>_fullflow.mat` the optical flow values stored as MATLAB variables (`U`, `V` for X and Y components of the vector field, and `Is` is the collapsed image stack with adjusted intensity.) 
* `<neuron_name>_fullflow_trace.avi` shows the optical flow overlaid by the neuron trace. 
* `<neuron_name>_ffd_color.avi` shows the stress in 2D collapsed neuron volume using the grid, where red means shrinkage, green means expansion, and yellow means no changes. 
* `<neuron_name>_fullflow_colorcode.avi` shows the stress in 2D collapsed neuron volume at every pixel using the color coded value (Red and yellow mean expansion in X and Y directions respectively, while green and purple show the shrinkage in X and Y directions respectively.)
* `<neuron_name>_regis.avi` shows the results of registration along the depth by the sequence of the collapsed registered neuron volume.
* `<neuron_name>_track.avi` shows the tracking results with the measurement of Calcium activity.
