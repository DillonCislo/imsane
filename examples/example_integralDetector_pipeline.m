%% ImSAnE Implementation of ectoderm extraction -- for Enrique Blanco
%
% We detect the midsurface of the folding ectoderm in the fly abdomen.
%
% This requires new tools compared to the original release of ImSAnE, 
% since (1) the surface cannot be described as a Monge form z(x,y) because
% of the 'multivalued' positions at the folds, and
% (2) 
% Additionally, while aspects of TubULAR have been incorporated into this
% workflow, the topology is disk-like, so we cannot simply rely on TubULAR.
% Instead, we use integralDetector (like in an ImSAnE analog of the
% TubULAR workflow) and then fit the surface with FoLD SurfER (Folding
% Surface & Lagrangian Deformation Extraction Resource, work in progress).
%
% This workflow was prepared by NPMitchell (npmitchell@kitp.ucsb.edu) and
% adapted from ImSAnE templates by Sebastian Streichan and Idse Heemskirk.

%% Set up ImSAnE paths
% Navigate to /path/to/imsane/setup.m
% Run setup.m

%% Initialize the project
%
% We start by creating an experiment object, which holds this metadata and 
% provides a frontend for a number of tasks such as loading the raw data.
% To construct the experiment object we need to pass dataDir, the directory 
% containing the raw data and projectDir, the directory where results of 
% the script will be saved.

clear all; close all;
cd /mnt/data/EnriqueBlanco/cropped_merged/

[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

dataDir = pwd;
projectDir = pwd;

xp = project.Experiment(projectDir, dataDir);

%%
% Next we set the metadata pertaining to the raw data files in the structure
% fileMeta. ImSAnE assumes that timepoints are saved individually and that 
% filenames in a timeseries are identical up to an integer specifying the 
% timepoint. Therefore we have
%
% * filenameFormat:               Filename, %u in the position of the integer
% * timePoints = [t1, t2, ..] :   List of times available. In this example 
% we have just a single time point 0.
% * stackResolution :             Stack resolution in micron.

fileMeta = struct();
fileMeta.dataDir = dataDir;
fileMeta.filenameFormat = 'merged_T%02d.tif' ; % 'ectoderm_c1_T%02d.tif';
fileMeta.timePoints = 0:19 ;
fileMeta.stackResolution = [0.69189 0.69189 1.0];   % resolution in um/pix 
fileMeta.rawStackSize = [210 182 50];               % #pix in X,Y, #frames in Z
fileMeta.nChannels = 1;

%% 
% In the structure expMeta we set general parameters for the surface
% analysis we will do. 
%
% * channelsUsed:   Which channels do we need.
% * channelColor:   Assign color to the channels, RGB = [1 2 3].
%                   In this example the first channel is E-cadherin, and 
%                   the second is actin. We want these in green and red,
%                   respectively.
% * dynamicSurface: Does the surface shape change with time?  
%                   For a single time point this is false. True is not yet
%                   supported.
% * jitterCorrection:   Not needed here.
% * detectorType:       Which type of surface detector will be used.
%                       We will look at only one side of the wing, which is
%                       a planar surface so we use planarDetector.
% * fitterType:         Which type of fitter will be used.
%                       We fit planar surfaces using Thin Plate Spline:
%                       tpsFitter.

expMeta = struct();
expMeta.description = 'Adult abdomen ectoderm';
expMeta.channelsUsed = [1 2];
expMeta.channelColor = [2 1];
expMeta.dynamicSurface = true;
expMeta.jitterCorrection = false;
expMeta.detectorType = 'surfaceDetection.integralDetector';
expMeta.fitterType = 'surfaceFitting.foldSurferFitter';

xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);

%%
% Finally we call initNew(), which reads the stack size from the first 
% available time point, then initializes fitter and detector and creates 
% fitOptions and detectOptions based on their defaults.

xp.initNew();

%% Load a time point from the data (OPTIONAL -- testing)
%
% Now that we have set up the project, we can load a time point.
% loadTime sets xp.currentTime, loads the stack into xp.stack 
% and resets detector and fitter with the default options for that time.
%
% We rescale to unit aspect ratio to make detection work better later and
% visualize in the right proportions.

% xp.loadTime(0);
% xp.rescaleStackToUnitAspect();

%% OPTIONAL -- testing
% xp.stack is not an array but a Stack object.
% The easy way to look at a slice through the data is using getSlice.

% imshow(xp.stack.getSlice('x', 300), []);

%% Detect the surface
%
% integralDetector.detectSurface detects the surface by minimizing a
% Chan-Vese free energy starting with an initial guess.
%

% Let's start out with the default options, but adjust some
detectOptions = xp.detector.defaultOptions;

% Set the diffusion constant for a Laplacian smoothing in MATLAB 
%   (WARNING: may require some dependencies from gptoolbox!)
detectOptions.smooth_with_matlab = 0.1 ;
detectOptions.tension = 1 ;

% We need pretty high resolution to resolve these folds, so let's subsample
% by only a small factor.
detectOptions.ssfactor = 2; 
detectOptions.enforceQuality = true ;

% Axis order of the data on disk
detectOptions.ilastikaxisorder = 'xyzc' ;

% Set some other params
detectOptions.include_boundary_faces = [0,0,0,0,1,0] ;  % we want a topological disk, not a sphere

%% Create subsampled h5 for training in iLastik
for tidx = [5, 10, 15, 20, 1:numel(xp.fileMeta.timePoints)]
    tp = xp.fileMeta.timePoints(tidx) ;
    fn = sprintf(xp.fileMeta.filenameFormat, tp) ;
    detectOptions.fileName = strrep(fn, '.tif', '.h5') ;
    detectOptions.timepoint = tp ;
    if ~exist(detectOptions.fileName, 'file')
        xp.loadTime(tp) ;
        xp.rescaleStackToUnitAspect() ;
        xp.setDetectOptions(detectOptions);
        xp.detector.prepareIlastik(xp.stack)
        % Note: 
        % If ssfactor ==4, (50/0.69189)/4 = 18, so we should expect 18 
        %   layers deep.
        %   And 973/4 = 243.25, so we should expect 243 layers wide/tall.
    else
        disp(['Found subsampled h5 on disk for t=' num2str(tp)])
    end
end

%% NOW TRAIN ON THE OUTPUTED h5 FILES IN ILASTIK TO FIND THE BASAL SURFACE
% OF THE TISSUE'S SEPTATE JUNCTIONS


%% Ensure that the output directory for meshes exists
if ~exist(detectOptions.meshDir, 'dir')
    disp(['Creating meshDir: ' detectOptions.meshDir])
    mkdir(detectOptions.meshDir)
end

%% Calling detectSurface runs the surface detector and creates the meshes
for tidx = 1:numel(xp.fileMeta.timePoints)
    tp = xp.fileMeta.timePoints(tidx) ;
    detectOptions.timepoint = tp ;
    detectOptions.niter = 200 ;
    detectOptions.niter0 = 2000 ;                                      
    detectOptions.niter = 100 ;     
    detectOptions.pre_pressure = -5 ;   
    detectOptions.post_pressure = 1 ;    
    detectOptions.tension = 0.01  ;     
    detectOptions.enforceQuality = false ;
    fn = sprintf(xp.fileMeta.filenameFormat, tp) ;
    detectOptions.fileName = strrep(fn, '.tif', '') ;
    xp.setDetectOptions(detectOptions);
    xp.detectSurface();
end

%% View all meshes (optional)
resolution = min(fileMeta.stackResolution) ;
% output directory for images as pngs
snapdir =fullfile(detectOptions.meshDir, 'snaps') ;
if ~exist(snapdir, 'dir')
    mkdir(snapdir)
end
for qq = 1:length(xp.fileMeta.timePoints)
    tp = xp.fileMeta.timePoints(qq) ;
    meshfn = fullfile(detectOptions.meshDir, sprintf('mesh_%06d.ply', tp)) ;
    mesh = read_ply_mod(meshfn) ;
    clf
    trisurf(triangulation(mesh.f, resolution*mesh.v), 'edgecolor', 'none')
    axis equal
    caxis([0, 20])
    zlim([0, 25]); 
    if qq == 1
        xlims = xlim ; ylims = ylim ;
    end
    xlim(xlims); ylim(ylims) ;
    grid off
    xlabel(['lateral position [' char(181) 'm]'])
    ylabel(['ap position [' char(181) 'm]'])
    zlabel(['dv position [' char(181) 'm]'])
    sgtitle(['t=' num2str(tp)])
    saveas(gcf, fullfile(snapdir, sprintf('mesh_%06d.png', tp)))
    pause(0.001)
end

%% (optional)
% We can also inspect a point cloud cross section over the data with
% detector.inspectQuality. In the pointCloud option, 'c' specifies the 
% color cyan.

xp.loadTime(19)
xp.stack.image.apply()
inspectOptions= struct('dimension', 'x', 'value', 32, 'pointCloud', 'c');
xp.detector.inspectQuality(inspectOptions, xp.stack);

%% (optional)
% Or we can look at the point cloud in 3d, with some subsampling factor.
ssfactor = 50;
xp.detector.pointCloud.inspect(ssfactor);

%% Parameterize the surface to a rectangular 2D domain minimizing Dirichlet energy
%
% By detecting the largest intensity jump along z for each x,y in the
% E-cad channel and filtering out local outliers we have found the apical
% surface of the disc proper cells. We can now fit a smooth surface
% representation to that.
%
% tpsFitter fits the pointcloud using a thin plate spline fit. It has the
% following options:
%
% * gridSize:     Size of grid on which to generate fitted surface
%               default [50 50], full size takes long.
% * smoothing:    TPS smoothing parameter (default 1000).

for tidx = 10:length(xp.fileMeta.timePoints)
    tp = xp.fileMeta.timePoints(tidx) ;
    disp(['Fitting surface for t=' num2str(tp)])
    
    % Skip this step if not inspecting the result    
    xp.loadTime(tp)
    IV = xp.stack.image.apply() ;
    
    % Now fit the surface (map to plane)
    fitOptions = struct('smoothing', 500, 'gridSize', [100 100]);
    meshFn = fullfile(detectOptions.meshDir, ...
        sprintf('mesh_%06d.ply', tp)) ;
    mesh = read_ply_mod(meshFn) ;
    xp.setFitOptions(fitOptions);
    %%
    options = struct() ;
    options.IV = IV ;
    try
    xp.fitSurface(mesh, options);

    % Save image to disk
    mkdir( fullfile(detectOptions.meshDir, 'pullbackImages'))
    imfn = fullfile(detectOptions.meshDir, 'pullbackImages', sprintf('max_%06d.tif', tp)) ;
    % imwrite(max(xp.fitter.fittedParam.im, [], 4), imfn)
    % imwrite(xp.fitter.fittedParam.im, imfn)
    % pause
    writeTiff5D(uint16((2^16-1)*xp.fitter.fittedParam.im), imfn, 16)
    catch
    end
end

%% NOW follow motion in pullback plane
error('have not written this part. Idea: run PIV in pullback to construct approx pathlines, then can track geometry of tissue as it deforms')

%% STILL NEED TO WORK ON THIS PART
if false 
    %%
    % We can visualize the result on a cross section with
    % fitter.inspectQuality.

    xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

    %%
    % The detector picks up the edge of the E-cad signal but the best read out
    % goes solidly through it so we want to move the surface down a little. For
    % this we use zEvolve, with a shift specified in pixels.

    shift = 12;
    xp.normallyEvolve(shift);
    xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

    %%
    % We now generate the Surface Of Interest. The charts to be generated are 
    % specified in xp.fitter.charts. In this case there is only one, called
    % 'uv'.
    xp.generateSOI();

    %% Pull back the data to the surface
    % 
    % We pull back the data to the SOI using pullbackStack.

    xp.SOI.pullbackStack(xp.stack, xp.currentROI, xp.currentTime);

    %%
    % To look at the pullback, we call the data field of the SOI at the right
    % time, and get a particular patch from that with getPatch. A patch is a 
    % part of a surface. In this case, there is only one called xy_index.
    % Then we get the data in some patch in a particular coordinate system with
    % getTransform. In this case there is only one coordinate system: xy.
    % What we get is an object not only holding the image data but also
    % metadata and methods to manipulate it. The actual data is obtained by
    % calling the method apply. This returns a cell array with entries for each
    % channel.

    % the first channel is Ecad
    channel = 1;

    discProperPatch = xp.SOI.data(tidx).getPatch('xy_index');
    discProperImage = discProperPatch.getTransform('xy').apply{channel};
    figure, imshow(discProperImage, [], 'InitialMagnification', 50);

    %% Save the result
    %
    % Finally we save the SOI using SOI.save. We set the following options:
    %
    % * dir:            The directory to save the SOI to.
    % * imwriteOptions: Pullbacks are saved to image files using imwrite, we
    % can pass options to change file format, compression etc. For example we
    % could change this option to
    % imwriteOptions = {'jp2', 'Mode', 'lossless'}; 
    % * make8bit:       Often absolute intensities don't matter and 8 bit offers
    % a large enough dynamic range. This options rescales the lookup table and
    % converts to 8 bit before saving.

    imwriteOptions = {'tif', 'Compression', 'deflate'};
    savedir = fullfile(scriptPath, 'discProperApicalSOI');

    options = struct(   'dir',              savedir,...
                        'imwriteOptions',   {imwriteOptions},...
                        'make8bit',         true);
    xp.SOI.save(options)

    %%
    % All metadata is saved in SOI.xml. Pullbacks, geometry and any other data
    % defined on the surface are saved to image files in subdirectories 
    % reflecting the structure of patches and coordinates explained earlier in 
    % this tutorial. We can reload a surface of interest with
    % SOI.load(directory)


end

