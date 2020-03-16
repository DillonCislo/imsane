classdef integralDetector_rawdata < surfaceDetection.surfaceDetector
    % Segmentation based on prediction maps from ilastik that maximize the
    % enclosed probability within a contiguous volume of pixels, given an
    % initial guess. A Chan-Vese energy functional is minimized, with three
    % terms: surface tension, pressure, and attachment energy which fixes
    % the boundary to high gradients in probabilities from ilastik
    % training.
    % 
    % properties: 
    %           defaultOptions : struct
    %           pointCloud : 
    %           options : struct
    % 
    % options:  channel : channel of the stack to use (RFP, etc)
    %           ssfactor : sub-sampling factor used in the ilastik
    %           classifier.
    %           fileName: Name of the ilastik prediction file. Typically
    %               ending with Probabilities.h5 in Ilastik v1.1
    %           zdim: cylinder axis in matlab coords, 2 = x
    
    %---------------------------------------------------------------------
    % license
    %---------------------------------------------------------------------

    % Copyright 2014 Idse Heemskerk, Sebastian Streichan, Noah Mitchell
    %
    % This file is part of ImSAnE.
    % 
    % ImSAnE is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    % 
    % ImSAnE is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    % 
    % You should have received a copy of the GNU General Public License
    % along with ImSAnE.  If not, see <http://www.gnu.org/licenses/>.
    %
    % We make an explicit exception giving permission to link this code
    % with GPL-incompatible matlab facilities.
    
    %---------------------------------------------------------------------
    % properties
    %---------------------------------------------------------------------
    
    properties (Constant)
        % default detector options
        defaultOptions = struct('channel', 1, ...
            'ssfactor', 4,... % subsampling factor: downsampling of raw data
            'niter', 40, ... % how many iterations before exit if no convergence
            'niter0', 90, ... % how many iterations before exit if no convergence for first timepoint
            'lambda1', 1, ...  % lambda1/lambda2 decides weight of inclusion/exclusion of interior/exterior
            'lambda2', 2, ...  % lambda1/lambda2 decides weight of inclusion/exclusion of interior/exterior
            'nu', 0.3, ... % float: how many pressure (dilation/erosion) steps per iteration
            'smoothing', 1,... % float: how many smoothing steps per iteration (can be <1)
            'post_nu', 3, ... % how many iterations to dilate (if positive) or erode (if negative) after convergence
            'post_smoothing', 2,... % how many iterations of smoothing after convergence
            'exit_thres', 1e-6, ... % convergence threshold: maximum difference between subsequent level sets upon which to exit algorithm ('close enough')
            'foreGroundChannel',2, ... % the index of the first dimension of the 4d input data (if 4d)
            'fileName',[], ... % the filename of h5 to train on
            'mslsDir', './msls_output/', ...  % the directory for all output data/images
            'ofn_ls', 'msls_apical_', ...  % the output filename for level sets
            'ofn_ply', 'mesh_apical_ms_', ... % the output filename for PLY files
            'ms_scriptDir', '/mnt/data/code/morphsnakes_wrapper/', ... % the directory containing run_morphsnakes.py
            'timepoint', 0, ... % which timepoint in the data to consider
            'zdim',2, ... % Which dimension is the z dimension
            'pre_nu', -5, ... % number of dilation/erosion passes for positive/negative values
            'pre_smoothing', 1, ... % number of smoothing passes before running MS
            'ofn_smoothply', 'mesh_apical_',... % the output file name (not including path directory)
            'mlxprogram', ... % the name of the mlx program to use to smooth the results
            './surface_rm_resample20k_reconstruct_LS3_1p2pc_ssfactor4.mlx',...
            'init_ls_fn', 'none', ... % the name of the initial level set to load, if any
            'run_full_dataset', false, ... % run MS on a time series, not just one file
            'radius_guess', -1, ... % radius of the initial guess sphere
            'dset_name', 'exported_data', ... % the name of the dataset to load from h5
            'clip', -1, ...  % if positive, value to clip the intensity of the raw data
            'save', false, ... % whether to save intermediate results
            'center_guess', 'empty_string', ... % xyz of the initial guess sphere ;
            'plot_mesh3d', false, ... % if save is true, plot intermediate results in 3d 
            'mask', 'none',... % filename for mask to apply before running MS
            'mesh_from_pointcloud', false, ... % Create a mesh from pointcloud
            'imsaneaxisorder', 'xyzc') % axis order relative to mesh axis order by which to process the point cloud prediction. To keep as mesh coords, use xyzc
            
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------    

    methods
         
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = integralDetector_rawdata()
            % Constructor
            %
            % radialEdgeDetector()
            
            this = this@surfaceDetection.surfaceDetector();
        end 
        
        % ------------------------------------------------------
        % surface detection
        % ------------------------------------------------------
        
        function detectSurface(this, ~)
            % Detect surface in the stack with preset options.
            %
            % detectSurface(ply_guess)
            %
            % Parameters
            % ----------
            % initial_guess : str path to h5 file of initial gues as binary implicit surface
            
            opts = this.options ;
            
            debugMsg(1, ['integralDetector.detectSurface() : channel='...
                num2str(opts.channel), ...
                ', ssfactor=' num2str(opts.ssfactor)...
                ', fileName =' num2str(opts.fileName)...
                ', zDim =' num2str(opts.zdim),'\n']);
            
            if isempty(opts.fileName)
                error('Please provide a regular prediction from ilastik in h5 format.');
            end
            % 
            % %---------------------------------
            % % Segmentation of a prediction map from ilastik. 
            % %---------------------------------
            % 
            % % load the exported data out of the ilastik prediction
            fileName = [opts.fileName,'.h5'];
            h5fileInfo = h5info(fileName);
            if strcmp(h5fileInfo.Datasets.Name,'inputData')
                im = h5read(fileName, '/inputData');
            elseif strcmp(h5fileInfo.Datasets.Name,'exported_data')
                im = h5read(fileName,'/exported_data');
            elseif strcmp(h5fileInfo.Datasets.Name,'volume')
                im = h5read(fileName,'/volume/prediction');
            else
                error(['Please provide a regular prediction from ilastik, either in', ...
                    'the format of version 1.1 or 0.5 (ie with exported_data as a dataset)']);
            end
            
            % size of prediction
            idxPerm = circshift(1:3, [1 -opts.zdim]);
            
            ySize = size(im, idxPerm(1));
            xSize = size(im, idxPerm(2));
            zSize = size(im, idxPerm(3));
            
            zmin = 1;
            zmax = size(im, idxPerm(3));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Identify the surface using the loaded probabilities here
            % Convert the current image to a level set using morphological snakes
            niter = opts.niter ;
            niter0 = opts.niter0 ;
            mslsDir = opts.mslsDir ;
            lambda1 = opts.lambda1 ;
            lambda2 = opts.lambda2 ;
            nu = opts.nu ;
            smoothing = opts.smoothing ;
            pre_nu = opts.pre_nu ;
            pre_smoothing = opts.pre_smoothing ;
            post_nu = opts.post_nu ;
            post_smoothing = opts.post_smoothing ;
            exit_thres = opts.exit_thres ;
            ofn_ply = opts.ofn_ply ; 
            ofn_ls = opts.ofn_ls ;
            channel = opts.channel ;
            ms_scriptDir = opts.ms_scriptDir ;
            % tpstamp = fileMeta.timePoints(first_tp);
            timepoint = opts.timepoint ;
            ofn_smoothply = opts.ofn_smoothply ;
            mlxprogram = opts.mlxprogram;
            init_ls_fn = opts.init_ls_fn ;
            clip = opts.clip ;
            % Run MS on a series of timepoints in a parent directory
            run_full_dataset = opts.run_full_dataset ;
            % Radius of initial guess if init_ls_fn does not exist or is
            % not supplied
            radius_guess = opts.radius_guess ;
            dset_name = opts.dset_name ;
            save = opts.save ;
            center_guess = opts.center_guess ;
            plot_mesh3d = opts.plot_mesh3d ;
            mask = opts.mask ;
            use_pointcloud = opts.mesh_from_pointcloud ;
            
            % Create the output dir if it doesn't exist
            if ~exist(mslsDir, 'dir')
                 mkdir(mslsDir)
            end
            
            % Convert the current image to a level set using morphological snakes
            %    '_Probabilities.h5']);
            
            scriptpath = fullfile(ms_scriptDir, 'run_morphsnakes.py') ;
            command = ['python ' scriptpath];
            % Check if we are running MS on a dataset or a single file
            if run_full_dataset 
                if ~strcmp(run_full_dataset, 'none') ...
                        && ~strcmp(run_full_dataset, '')
                    use_dataset_command = true ;
                else
                    use_dataset_command = false ;
                end
            else 
                use_dataset_command = false ;
            end
            
            if use_dataset_command
                % User has elected to run as a dataset, so pass a directory
                % with _Probabilities.h5 files to run on.
                command = [command ' -dataset'] ;
                command = [command ' -i ' run_full_dataset ] ;
                msls_mesh_outfn = ofn_ply;
                ls_outfn = ofn_ls;
            else
                % We are running MS on a single file. Give the filename of
                % the ilatik output run on filename.h5
                prob_infn = [opts.fileName, '.h5'] ;
                command = [command ' -i ' prob_infn ] ;
                msls_mesh_outfn = [ofn_ply, num2str(timepoint, '%06d'), '.ply'];
                ls_outfn = [ofn_ls, num2str(timepoint, '%06d'), '.h5'];
            end
            outputLs = fullfile(mslsDir, ls_outfn) ;
            outputMesh = fullfile(mslsDir, msls_mesh_outfn) ;
            command = [command ' -o ' mslsDir ] ;
            command = [command ' -prenu ' num2str(pre_nu) ' -presmooth ' num2str(pre_smoothing)] ;
            command = [command ' -ofn_ply ' msls_mesh_outfn ' -ofn_ls ' ls_outfn];
            command = [command ' -l1 ' num2str(lambda1) ' -l2 ' num2str(lambda2) ] ;
            command = [command ' -nu ' num2str(nu) ' -postnu ' num2str(post_nu) ];
            command = [command ' -channel ' num2str(channel) ] ;
            command = [command ' -smooth ' num2str(smoothing) ];
            command = [command ' -postsmooth ' num2str(post_smoothing) ];
            command = [command ' -exit ' num2str(exit_thres, '%0.9f') ];
            command = [command ' -dset_name ' dset_name ] ;
            command = [command ' -dtype h5 -clip ' num2str(clip) ] ;
            command = [command ' -permute xyz' ] ;
            command = [command ' -ss ' num2str(opts.ssfactor)] ;
            if plot_mesh3d
                command = [command ' -plot_mesh3d' ] ;
            end
            if save
                command = [command ' -save'] ;
            end
            if ~strcmp(center_guess, 'empty_string')
                command = [command ' -center_guess ' center_guess ];
            end
            
            % get mask if supplied
            if ~strcmp(mask, 'none') && ~strcmp(mask, 'empty_string')
                command = [ command ' -mask ' mask ] ;
            end
            
            
            if radius_guess > 0
                command = [command ' -rad0 ' num2str(radius_guess)] ;
            end
            
            % Check if previous time point's level set exists to use as a seed
            % First look for supplied fn from detectOptions. 
            % If not supplied (ie init_ls_fn is none or empty string, then
            % seek previous timepoint output from MS algorithm.
            if strcmp(init_ls_fn, 'none') || strcmp(init_ls_fn, '')
                % User has NOT supplied fn from detectOptions
                init_ls_fn = [mslsDir 'msls_apical_', ...
                    num2str(timepoint - 1, '%06d' ) '.h5'] ;
            end
            
            if exist(fullfile(mslsDir, init_ls_fn), 'file')
                % It does exist. Use it as a seed (initial level set)
                command = [command ' -init_ls ', ...
                    fullfile(mslsDir, init_ls_fn), ...
                    ' -n ' num2str(niter) ] ;
            else
                % The guess for the initial levelset does NOT exist, so use
                % a sphere for the guess. 
                disp('Using default sphere for init_ls')
                command = [command ' -n ' num2str(niter0)];
            end
            disp(['prepared command = ', command])
            if ~exist(outputMesh, 'file')
                % Either copy the command to the clipboard
                clipboard('copy', command);
                % or else run it on the system
                system(command) 
            else
                disp(['output PLY already exists: ', msls_mesh_outfn])
            end

             % Clean up mesh file for this timepoint using MeshLab --------
            msls_mesh_outfn = [ofn_ply, num2str(timepoint, '%06d' ), '.ply'];
            PCfile = fullfile( mslsDir, msls_mesh_outfn );
            mesh_outfn = [ofn_smoothply, ...
                num2str(timepoint, '%06d'), '.ply'];
            outputMesh = fullfile(mslsDir, mesh_outfn);

            if ~exist( outputMesh, 'file')
                command = ['meshlabserver -i ' PCfile ' -o ' outputMesh, ...
                           ' -s ' mlxprogram ' -om vn'];
                if exist(PCfile, 'file')
                    % Either copy the command to the clipboard
                    clipboard('copy', command);
                    % or else run it on the system
                    disp(['running ' command])
                    system(command) 
                else
                    error(['unsmoothed PLY not found: ' PCfile])
                end
            else
                disp(['t=', num2str(timepoint) ': smoothed mesh file found, loading...'])
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~exist(outputMesh, 'file')
                error(['Meshlab did not create the file: ' outputMesh])
            end
            
            disp(['reading PLY ', outputMesh])
            tmp = read_ply_mod(outputMesh);
            vv = tmp.v ;
            points = struct('x', vv(:, 1), 'y', vv(:, 2), 'z', vv(:, 3));
            
            % don't forget to rescale the data;
            x = cat(1,points.x);
            y = cat(1,points.y);
            z = cat(1,points.z);     
            pointCloud = [x,y,z];            
            
            %--------------------------------------------------
            % scale point cloud to full size and set alignment
            %--------------------------------------------------
            
            largePC = zeros(size(pointCloud));
            for i = 1:3
                largePC(:,i) = (pointCloud(:,i)-1)*opts.ssfactor + 1;
            end
            
            largePC = sortrows(largePC, 3);
            
            debugMsg(2, 'setting ROI alignment based on zDim\n');
            
            ROI = surfaceDetection.RegionOfInterest(eye(4), eye(4));
            if opts.zdim == 2 % i.e zp == x bc x is the second index
                ROI.setAxes([0 1 0], [0 0 1], [1 0 0]);
            elseif opts.zdim == 1    
                ROI.setAxes([0 0 1], [1 0 0], [0 1 0]);
            end
            % we also want to set the ranges of the ROI, i.e. some bounding 
            % box that contains the pointcloud which will be used by the
            % fitter to determine its initial domain
            xpRange = [1, xSize*opts.ssfactor];
            ypRange = [1, ySize*opts.ssfactor];
            zpRange = [1, zSize*opts.ssfactor];
            ROI.setRanges(xpRange,ypRange,zpRange);
            this.pointCloud = surfaceDetection.PointCloud(largePC,ROI);
            this.pointCloud.determineROI(15);

        end
        
        % ------------------------------------------------------
        % prepare data for ilastik segmentation
        % ------------------------------------------------------
        
        function prepareDownsampled(this,stack)
            % Prepaere stack for Ilastik segmetnation. This outputs an h5
            % file of subsampled intensity data on which to train.
            %
            % prepareIlastik(stack)
            
            % Accoring to the specified options sub-sample the stack and
            % save it for analysis with ilastik. 
            
            opts = this.options;
            im = stack.image.apply();
            fileName = [opts.fileName,'.h5'];
            
            if exist(fileName,'file')
                disp(['deleting downsampled data on disk: ' fileName])
                delete(fileName)
            end    
                  
            dsetName = '/inputData';
            
            for c = 1 : length(im)                
                image(:,:,:,c) = im{c}(1:opts.ssfactor:end,1:opts.ssfactor:end,1:opts.ssfactor:end);
            end
            
            disp(['Writing downsampled dataset to disk: ' fileName])
            if ndims(image)==4
                h5create(fileName,dsetName,size(image));
                h5write(fileName,dsetName,image);
            else
                h5create(fileName,dsetName,size(image));
                h5write(fileName,dsetName,image);
            end
            
        end
        
        % ------------------------------------------------------
        % Check point cloud against original image
        % ------------------------------------------------------
        function inspectQuality(this, inspectOpts, stack)
            %   inspect quality of fit in single slice in dimension specified 
            %   by options and display image.
            %
            %   inspectQuality(inspectOpts, stack)
            %
            %   override inspectQuality to deal with subsampled point cloud
            
            ssfac = this.options.ssfactor;
            
            if ssfac ~= 1 && rem(inspectOpts.value-1, ssfac) ~= 0
                inspectOpts.value = round(inspectOpts.value/ssfac)*ssfac + 1;
                disp([]);
                debugMsg(2, ['WARNING: sub-sampled point cloud taken from different'...
                    ' plane, may look a little off\n']);
            end

            inspectQuality@surfaceDetection.surfaceDetector(this, inspectOpts, stack);
        end
        
    end
end
