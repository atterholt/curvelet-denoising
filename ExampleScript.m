% Script to get people acquainted with curvelet denoising
% James Atterholt Caltech 2021

% If you use this work, please use the reference:
% Atterholt, J., Zhan, Z., Shen, Z., Li, Z., (Under Review). A unified 
% wavefield-partitioning approach for distributed acoustic sensing.

clear

% Add path to the CurveLab Toolbox (put your own path in here). You can
% download this software package at www.curvelet.org
addpath("~/SeismoPrograms/CurveLab-2.1.3/fdct_wrapping_matlab/")

% Add path to your favorite seismic colormap. Here I use bluewhitered.
% Link: https://www.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered
addpath("~/SeismoPrograms/ColorMaps/")

% Table of contents:
% 1: Loading in the data
% 2: Performing the preprocessing
% 3: Understanding the coherent filtering
% 4: Understanding the stochastic filtering
% 5: Applying the curvelet filter
% 6: Visualizing your results

% Figures:
% 1: Raw data
% 2: Preprocessed raw data
% 3: Table of velocity bounds from which to choose
% 4: Raw noise vector
% 5: Preprocessed noise vector
% 6: Coherent filtered data
% 7: Stochastic filtered data
% 8: Both coherent and stochastic filtered data

% Things you need to define (if you are modifying this to your own data):
% 1: dt (sampling rate) - section 1
% 2: ds (station spacing) - section 1
% 3: BadSta (bad stations to be removed as preprocessing) - section 1
% 4: nbangles (number of angles [polar tiling parameter]) - section 3
% 5: nbscales (number of scales [polar tiling parameter]) - section 3

% References:
% E. J. Candes, L. Demanet, D. L. Donoho, L. Ying, Fast Discrete Curvelet 
%   Transforms, 2005.
% E. J. Candes, D. L. Donoho, New Tight Frames of Curvelets and Optimal 
%   Representations of Objects with Smooth Singularities, 2002.

%% Loading in the data

% Load in the example data as EQ_raw
load("EQ_raw.mat")

% Define a priori variables
dt = (1/250); % s (sampling rate)
ds = 8; % m (station spacing)
f = 1/dt; % hz (frequency)
k = 1/ds; % 1/m (spatial frequency)

% Get distance and time
t = (1:size(EQ_raw,1))*dt;
d = (1:size(EQ_raw,2))*ds;

% Define bad and good stations for this array (hand picked)
BadSta = [1:24 83:90 143:147 256:271 455:473 663:670 843:872 1060:1067];
GoodSta = setdiff(1:size(EQ_raw,2),BadSta);

% Inspect the data
figure(1)
imagesc(d,t,EQ_raw);
caxis([-1 1]*10);
colormap(bluewhitered(256))
xlabel("Distance (m)")
ylabel("Time (s)")
title("Raw Data")

%% Performing the preprocessing

% Get rid of the bad stations
EQ_raw_b = EQ_raw(:,GoodSta);

% Apply a median filter (a,b,c selected empirically, see function for details)
EQ_raw_c = MedianFilter(EQ_raw_b,50,5,10);

% Get new spatial axis (you've now gotten rid of a bunch of stations)
d2 = (1:size(EQ_raw_c,2))*ds;

% Visualize preprocessed data
figure(2)
imagesc(d2,t,EQ_raw_c);
caxis([-1 1]*10);
colormap(bluewhitered(256))
xlabel("Distance (m)")
ylabel("Time (s)")
title("Raw Data (Preprocessed)")

%% Understanding the coherent filtering

% Define the number of angles at the coursest scale. I've chosen to start
% coherent filtering one scale above the coursest scale, because the
% coursest scale is really too large to deal with any seismological
% information. If it ain't broke, don't fix it, eh? This has implications
% for the fineness of wedges you can choose to filter. If you want, say, 16
% angles at the coursest scale, then you'd start coherent filtering at the
% scale above, which will always have 2x as many angles, in this case 32
% angles. This gets confusing because half of the wedges are mirrors of the
% other half (for our purposes), so that makes for just 16 velocity ranges.
% This really sounds like a mess, but all you're doing is defining how
% finely you want to tile the polar tiling of the FK frame. nbangles is
% good for most DAS arrays at 16, and nbscales (the number of scales, or
% number of bins away from the origin in FK) is generally good at 8. Note
% that your choice for nbangles must be a factor of 4.
nbangles = 16;
nbscales = 8;

% Assign velocity ranges to wedges (these are assigned according to
% velocities defined by the wedges in FK space). Each wedge has the same 
% base-length and the number of such wedges are assigned according to 
% nbangles (see above). If you want more precision in your velocity 
% filtering, you can always add more wedges (you just need to have wedges
% with factors of 4). This script should work for an arbitrary nbangles
% (given that your nbangles is a factor of 4), but FYI you can understand
% how these are defined by breaking up your FK space into 2*nbangles
% triangles where each triangle has the same base width (the base being the
% side on the edge of the FK plane). Then the other two sides of the
% triangle define slopes in FK space, which have equivalent velocities. So,
% what we're computing below are slopes! 
v_bounds = zeros((nbangles/2)+1,1);
halfway = nbangles/4;
v_bounds(halfway+1,1) = (f)/(k);
for i = 1:(halfway)
    v_bounds(i,1) = (((i-1)/halfway)*f)/(k);
    v_bounds(halfway+i+1,1) = (f)/(((halfway-i)/halfway)*k);
end

% There are integers assigned to each wedge AND each direction (west-going 
% waves and east-going waves). The starting point of these integers is kind 
% of arbitrary so I'll just define the integers for you. Just note that
% these wedges' bounds start from the top left corner and rotate clockwise 
% assuming you've plotted your FK spaces oriented as (x,y) = (k,f). So,
% each wedge gets its own integer. Then you can reference the velocity
% range to zero it with that integer.
ints = 1:nbangles; 
dirs = strings(nbangles,1);
for i = 1:nbangles
   if (i <= halfway) || (i > (3*halfway))
       dirs(i) = "W";
   else
       dirs(i) = "E";
   end
end
v_lows = [(halfway+1):halfway*2 halfway*2:-1:(halfway+1) halfway:-1:1 1:halfway];
v_highs = v_lows + 1;

% Now define a string array to describe integer choices. This just puts the
% velocity bounds with their assigned integers
int_choices = strings(nbangles,1);
for i = 1:size(int_choices,1)
    int_choices(i) = strcat(num2str(ints(i)),": ",num2str(round(v_bounds(v_lows(i))))," to ",num2str(round(v_bounds(v_highs(i))))," m/s ",dirs(i));
end

% Plot these for visualization purposes (ignore this, it's gross), but the
% result is a table with the integer associated with each velocity range.
% This is how you'll choose the velocity ranges you'd like to filter out.
% You put the integer associated with the velocity range in the vector
% vels_2_zero below.
plot_locs = zeros(nbangles,2);
irange = 4;
jrange = nbangles/4;
for i = 1:irange
    for j = 1:jrange
        n = (i-1)*jrange + j;
        plot_locs(n,1) = j;
        plot_locs(n,2) = -i+5;
    end
end
wedges = figure(3);
wedges.Position(3:4) = [650*(jrange/4) 200];
text(plot_locs(:,1),plot_locs(:,2),int_choices)
xlim([.5 halfway+1])
ylim([0 4+1])

% Okay so, say we want to get rid of all velocities under 1000 m/s to 
% remove the traffic noise, and if we've set our nbangles to 16, then by 
% our chart, we want to zero wedges 11, 12, 13, and 14:
vels_2_zero = [11,12,13,14];

% To do the same filtering with nbangles = 32, then uncomment this:
% vels_2_zero = [21,22,23,24,25,26,27,28];

% To keep our results real on our inverse transform, we need our zeroing to
% be symmetric. We accomplish this by taking advantage of the fact that we
% can find the mirror wedge by adding a constant value, which happens to be
% the nbangles.
vels_2_zero = [vels_2_zero , (vels_2_zero + nbangles)];

%% Understanding the stochastic noise filtering

% Stochastic noise filtering can be done empirically in a couple of ways.
% If you have a pure noise window, you can establish thresholds using that
% window. We'll do that first

% Load in the noise window and process it as you did the raw data
load("Noise_Matrix.mat")

% Inspect the noise matrix
figure(4)
imagesc(d,t,Noise_Matrix);
caxis([-1 1]*10);
colormap(bluewhitered(256))
xlabel("Distance (m)")
ylabel("Time (s)")
title("Noise Matrix")

% Get rid of the bad stations
Noise_Matrix_b = Noise_Matrix(:,GoodSta);

% Apply a median filter
Noise_Matrix_c = MedianFilter(Noise_Matrix_b,50,5,10);

% Visualize the preprocessed noise (that sounds stupid, doesn't it?)
figure(5)
imagesc(d2,t,Noise_Matrix_c);
caxis([-1 1]*10);
colormap(bluewhitered(256))
xlabel("Distance (m)")
ylabel("Time (s)")
title("Noise Matrix (Preprocessed)")

% Now, transform this over to the curvelet domain. Because this is the
% first curvelet transform we've done, I'll just explain really quickly
% what each argument of this function does:
% 1st: The matrix that you are transforming
% 2nd: Real or complex valued curvelets (0 = complex, 1 = real)
% 3rd: Objects at the finest scale, I choose wavelets for expediency (1 =
% curvelets, 2 = wavelets)
% 4th and 5th: These define the structure of the polar tiling you are
% using (4th: number of scales, 5th: number of angles at the coursest scale). 
% We end up in the weeds here. I choose 8 and 16 because those are
% optimized sizes for the example and most DAS arrays. Number of angles may
% be increased for more finely compartmentalized velocity bounds.
% Note that I've chosen the 2nd and 3rd arguments for you. This is beyond
% the scope of this script. If you'd like more details about this I'd
% direct you to the original paper on the fast discrete curvelet transform
C = fdct_wrapping(Noise_Matrix_c,1,2,nbscales,nbangles);

% Now, since we're going to perform soft thresholding on these curvelets,
% we need to define thresholds. We can define the thresholds using our
% noise window by using the curvelet coefficients of the noise window. We
% get these using the following code:
% Compute norm of curvelets (exact)
E_noise_vector = cell(size(C));
for s=1:(length(C))
    E_noise_vector{s} = cell(size(C{s}));
    for w=1:length(C{s})
    
        % Allocate for a vector of all coefficients in a wedge
        WedgeCoeffs = numel(C{s}{w}(:));
                
        % Vectorize all coefficients in a scalar threshold
        WedgeLinear = abs(reshape(C{s}{w},WedgeCoeffs,1));
        
        % Get the 95% of the vector (can be more or less conservative here)
        E_noise_vector{s}{w} = prctile(WedgeLinear,95);
    
    end
    
end

% Another option if you don't have a noise vector is to perform a knee
% point based thresholding (curvelets with information will have high
% amplitudes and curvelets without information will have low amplitudes, so
% you can take the iECDF and notice that there are a few curvelets with
% high amplitudes. The knee point of this plot gives you an estimate of
% where the transition between non-information-having and
% information-having curvelets is). This is old, because I used the noise
% vector approach in the paper, so there may be better approaches here, but
% this is one option for you to use.

% Again, go to the curvelet domain which is captured in C
C = fdct_wrapping(EQ_raw_c,1,2,nbscales,nbangles);

E_knee_points = cell(size(C));

% Run loop to determine threshold values
for s=1:(length(C))
    E_knee_points{s} = cell(size(C{s}));
    for w=1:length(C{s})
        
        % Allocate for a vector of all coefficients in a wedge
        WedgeCoeffs = numel(C{s}{w}(:));
    
        % Vectorize all coefficients in a scalar threshold
        WedgeLinear = abs(reshape(C{s}{w},WedgeCoeffs,1));
        
        % Take the empirical cumulative dist of your new vector
        [F,x] = ecdf(abs(WedgeLinear));
    
        % Then compute the slope from the first to last element in the 
        % inverse ECDF
        m = ((x(end)-x(1))/(F(end)-F(1))).*F;
    
        % Use the slope to find the knee point
        tiltedplot = x-m;
        idx = (tiltedplot == min(tiltedplot));
    
        % Now set the knee point as the threshold. Weaken this threshold by
        % some factor (here it's 0.25 to be conservative). The
        % multiplication factor is a dial for you to test to see how
        % conservative or unconservative you wish to be. 
        E_knee_points{s}{w} = x(idx)*0.25;
        
    end
end

%% Applying the curvelet filter

% Here we'll apply our curvelet filtering approach with only coherent noise
% filtering, only stochastic noise filtering, and both. Feel free to pick
% your poison. This is all baked into a function called CurveletDenoising.
% Note you will need to edit the coherent denoising part of
% CurveletDenoising if you want to alter the polar tiling structure. It is
% pretty clear what I do in that part, that is, I propagate the zeroing of
% wedges up the scales so that I end up zeroing one big wedge. The
% nonisotropic scaling adds some complexity to this algorithm.

% Argument 1: Matrix you want filtered
% Argument 2: Stochastic noise filter cell (shape of intended curvelet transform)
% Argument 3: Coherent noise filter vector (vector w/ number of elements
% equal to the number of angles at the 2nd coarsest scale)
% Argument 4: Choice = 0 (just stochastic) or 1 (just coherent) or 2 (both)
% Note: If choice = 0 or 1 just put in an empty vector for either E or v
% Argument 5: The number of scales used (see definition above)
% Argument 6: The number of angles used (see definition above)

% Do just the coherent noise filtering
EQ_coherent = CurveletDenoising(EQ_raw_c,E_noise_vector,vels_2_zero,1,nbscales,nbangles);

% Do just the stochastic noise filtering. Here I use the empirical noice
% vector E. But if you'd like to try the knee point method, just replace
% E_noise_vector with E_knee_points
EQ_stochastic = CurveletDenoising(EQ_raw_c,E_noise_vector,vels_2_zero,0,nbscales,nbangles);

% Do both
EQ_both = CurveletDenoising(EQ_raw_c,E_noise_vector,vels_2_zero,2,nbscales,nbangles);

%% Visualizing your results

% Plot just the coherent filtered data
figure(6)
imagesc(d2,t,EQ_coherent);
caxis([-1 1]*10);
colormap(bluewhitered(256))
xlabel("Distance (m)")
ylabel("Time (s)")
title("Coherent Filtered")

% Plot just the stochastic filtered data
figure(7)
imagesc(d2,t,EQ_stochastic);
caxis([-1 1]*10);
colormap(bluewhitered(256))
xlabel("Distance (m)")
ylabel("Time (s)")
title("Stochastic Filtered")

% Plot the data filtered for both coherent and stochastic noise
figure(8)
imagesc(d2,t,EQ_both);
caxis([-1 1]*10);
colormap(bluewhitered(256))
xlabel("Distance (m)")
ylabel("Time (s)")
title("Both Filtered")









