classdef integralDetector < surfaceDetection.surfaceDetector
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
    %           foreGroundChannel: the labeled used in ilastik as
    %               foreground default is 2. (1 would then
    %               be the first channel of the ilastik training).
    %               Note that this is 1-indexed, so 1 would be the first
    %               channel.
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
            'mlxprogram', ... % the name of the mlx program to use to smooth the results. Note that if mesh_from_pointcloud==true, should take obj as input and mesh as output.
            './surface_rm_resample20k_reconstruct_LS3_1p2pc_ssfactor4.mlx',...
            'init_ls_fn', 'none', ... % the name of the initial level set to load, if any
            'run_full_dataset', false, ... % run MS on a time series, not just one file
            'radius_guess', -1, ... % radius of the initial guess sphere
            'dset_name', 'exported_data', ... % the name of the dataset to load from h5
            'save', false, ... % whether to save intermediate results
            'center_guess', 'empty_string', ... % xyz of the initial guess sphere ;
            'plot_mesh3d', false, ...  % if save is true, plot intermediate results in 3d 
            'dtype', 'h5', ... % h5 or npy: use hdf5 or numpy file format for input and output ls
            'mask', 'none', ... % filename for mask to apply before running MS
            'mesh_from_pointcloud', false) ; % use a pointcloud from the marching cubes algorithm rather than a mesh to create smoothed mesh
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = integralDetector()
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
            % initial_guess : str path to npy file of initial gues as binary implicit surface
            
            opts = this.options ;
            
            debugMsg(1, ['integralDetector.detectSurface() : channel='...
                num2str(opts.channel), ...
                ', ssfactor=' num2str(opts.ssfactor)...
                ', fileName =' num2str(opts.fileName)...
                ', foreGroundChannel =' num2str(opts.foreGroundChannel)...
                ', zDim =' num2str(opts.zdim),'\n']);
            
            if isempty(opts.fileName)
                error('Please provide a regular prediction from ilastik in h5 format.');
            end
            
            %---------------------------------
            % Segmentation of a prediction map from ilastik.
            %---------------------------------
            
            % load the exported data out of the ilastik prediction
            fileName = [opts.fileName,'_Probabilities.h5'];
            disp(['Reading h5 file: ' fileName])
            h5fileInfo = h5info(fileName);
            if strcmp(h5fileInfo.Datasets.Name,'exported_data')
                file = h5read(fileName,'/exported_data');
            elseif strcmp(h5fileInfo.Datasets.Name,'volume')
                file = h5read(fileName,'/volume/prediction');
            else
                error(['Please provide a regular prediction from ilastik, either in', ...
                    'the format of version 1.1 or 0.5 (ie with exported_data as a dataset)']);
            end
            
            foreGround = opts.foreGroundChannel;
            
            % ilastik internally swaps axes. 1: class, 2: y, 3: x 4 : z
            pred = permute(file,[3,2,4,1]);
            pred = pred(:,:,:,foreGround);
            pred = uint8(255*pred);
            
            % size of prediction
            idxPerm = circshift(1:3, [1 -opts.zdim]);
            
            ySize = size(pred, idxPerm(1));
            xSize = size(pred, idxPerm(2));
            zSize = size(pred, idxPerm(3));
            
            zmin = 1;
            zmax = size(pred, idxPerm(3));
            
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
            % Run MS on a series of timepoints in a parent directory
            run_full_dataset = opts.run_full_dataset ;
            % Radius of initial guess if init_ls_fn does not exist or is
            % not supplied
            radius_guess = opts.radius_guess ;
            save = opts.save ;
            center_guess = opts.center_guess ;
            plot_mesh3d = opts.plot_mesh3d ;
            dtype = opts.dtype ; 
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
                prob_infn = [opts.fileName, '_Probabilities.h5'] ;
                command = [command ' -i ' prob_infn ] ;
                msls_mesh_outfn = [ofn_ply, num2str(timepoint, '%06d'), '.ply'];
                ls_outfn = [ofn_ls, num2str(timepoint, '%06d'), '.', dtype];
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
            command = [command ' -channel ' num2str(foreGround - 1) ] ;
            command = [command ' -dtype ' dtype ] ; 
            if plot_mesh3d
                command = [command ' -plot_mesh3d' ] ;
            end
           
            if save
                command = [command ' -save'] ;
            end
            
            if ~strcmp(center_guess, 'empty_string')
                command = [command ' -center_guess ' center_guess ];
            end
            
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
            
            disp(init_ls_fn)
            if strcmp(init_ls_fn, 'none') || strcmp(init_ls_fn, '')
                % User has NOT supplied fn from detectOptions
                init_ls_fn = ['msls_apical_', ...
                    num2str(timepoint - 1, '%06d' ) '.' dtype] ;
            end
            
            disp(init_ls_fn)
            if exist(fullfile(mslsDir, init_ls_fn), 'file')
                % It does exist. Use it as a seed (initial level set)
                disp('running using initial level set')
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

            % Here use the boundary mesh from marching cubes to make a
            % smooth mesh

            if use_dataset_command
                % find all ms_...ply files in mslsDir, and smooth them all
                files_to_smooth = dir(fullfile(mslsDir, [ofn_ply '*.ply'])) ;
                lsfns_to_smooth = dir(fullfile(mslsDir, [ls_outfn '*' dtype])) ;
                
                if length(lsfns_to_smooth) ~= length(files_to_smooth)
                    error('The number of output levelsets does not equal the number of output unsmoothed PLYs. These must match.')
                end
                for i=1:length(files_to_smooth)
                    msls_mesh_outfn = files_to_smooth(i).name ;
                    PCfile = fullfile( mslsDir, msls_mesh_outfn );
                    % Note that LS file is outputLs ;
                    split_fn = strsplit(msls_mesh_outfn, ofn_ply) ;
                    extension_outfn = split_fn{2} ;
                    base_outfn_for_pointcloud = ofn_ply ;
                    mesh_outfn = [ofn_smoothply, extension_outfn] ;
                    outputMesh = fullfile(mslsDir, mesh_outfn);

                    if ~exist( outputMesh, 'file')
                        if use_pointcloud
                            % Use the pointcloud from the level set rather than the
                            % boundary mesh from marching cubes
                            %----------------------------------------------------------------------
                            % Extract the implicit level set as a 3D binary array
                            %----------------------------------------------------------------------

                            % The file name of the current time point
                            ls_outfn_ii = lsfns_to_smooth(i).name ;
                            % The 3D binay array
                            bwLS = h5read( ls_outfn_ii, '/implicit_levelset' );

                            % Extract the (x,y,z)-locations of the level set boundary (in pixel
                            % space)
                            bwBdyIDx = bwperim( bwLS );

                            clear bwBdy
                            [ bwBdy(:,1), bwBdy(:,2), bwBdy(:,3) ] = ind2sub( size(bwLS), ...
                                find(bwBdyIDx) );

                            %----------------------------------------------------------------------
                            % Create output mesh
                            %----------------------------------------------------------------------

                            % Write the points to a .obj file as a point cloud for ouput to Meshlab
                            clear OBJ
                            OBJ.vertices = bwBdy;
                            OBJ.objects(1).type='f';
                            OBJ.objects(1).data.vertices=[];
                            
                            pointcloud_fn = [base_outfn_for_pointcloud '.obj'] ;
                            disp(['Writing point cloud ' pointcloud_fn]);
                            write_wobj(OBJ, pointcloud_fn );

                            % Run the meshlab script
                            system( ['meshlabserver -i ' pointCloudFileName, ...
                                ' -o ' outputMesh, ...
                                ' -s ' mlxprogram ' -om vn']);
                        else
                            if use_pointcloud
                                % Use the pointcloud from the level set rather than the
                                % boundary mesh from marching cubes
                                %----------------------------------------------------------------------
                                % Extract the implicit level set as a 3D binary array
                                %----------------------------------------------------------------------

                                % The file name of the current time point is ls_outfn 
                                % The 3D binay array
                                bwLS = h5read( outputLs, '/implicit_levelset' );

                                % Extract the (x,y,z)-locations of the level set boundary (in pixel
                                % space)
                                bwBdyIDx = bwperim( bwLS );

                                clear bwBdy
                                [ bwBdy(:,1), bwBdy(:,2), bwBdy(:,3) ] = ind2sub( size(bwLS), ...
                                    find(bwBdyIDx) );

                                %----------------------------------------------------------------------
                                % Create output mesh
                                %----------------------------------------------------------------------

                                % Write the points to a .obj file as a point cloud for ouput to Meshlab
                                clear OBJ
                                OBJ.vertices = bwBdy;
                                OBJ.objects(1).type='f';
                                OBJ.objects(1).data.vertices=[];

                                pointcloud_fn = [base_outfn_for_pointcloud '.obj'] ;
                                disp(['Writing point cloud ' pointcloud_fn]);
                                write_wobj(OBJ, pointcloud_fn );

                                % Run the meshlab script
                                system( ['meshlabserver -i ' pointCloudFileName, ...
                                    ' -o ' outputMesh, ...
                                    ' -s ' mlxprogram ' -om vn']);
                            else
                                % Use the marching cubes mesh surface to smooth
                                command = ['meshlabserver -i ' PCfile ' -o ' outputMesh, ...
                                    ' -s ' mlxprogram ' -om vn'];
                                % Either copy the command to the clipboard
                                clipboard('copy', command);
                                % or else run it on the system
                                disp(['running ' command])
                                system(command)
                            end
                        end
                    else
                        disp(['t=', num2str(timepoint) ': smoothed mesh file found...'])
                    end

                end
            else
                msls_mesh_outfn = [ofn_ply, num2str(timepoint, '%06d' ), '.ply'];
                PCfile = fullfile( mslsDir, msls_mesh_outfn );
                mesh_outfn = [ofn_smoothply, num2str(timepoint, '%06d'), '.ply'];
                outputMesh = fullfile(mslsDir, mesh_outfn);

                if ~exist( outputMesh, 'file')
                    command = ['meshlabserver -i ' PCfile ' -o ' outputMesh, ...
                        ' -s ' mlxprogram ' -om vn'];
                    % Either copy the command to the clipboard
                    clipboard('copy', command);
                    % or else run it on the system
                    disp(['running ' command])
                    system(command)
                else
                    disp(['t=', num2str(timepoint) ': smoothed mesh file found, loading...'])
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if exist(outputMesh, 'file')
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
            else
                error(['Output mesh from detector not found! Sought: ' outputMesh])
            end
            
        end
        
        % ------------------------------------------------------
        % prepare data for ilastik segmentation
        % ------------------------------------------------------
        
        function prepareIlastik(this,stack)
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
                delete(fileName)
            end
            
            dsetName = '/inputData';
            
            for c = 1 : length(im)
                image(:,:,:,c) = im{c}(1:opts.ssfactor:end,1:opts.ssfactor:end,1:opts.ssfactor:end);
            end
            if ndims(image)==4
                h5create(fileName,dsetName,[size(image,2) size(image,1), size(image,3) size(image,4)]);
                h5write(fileName,dsetName,permute(image,[2 1 3 4]));
            else
                h5create(fileName,dsetName,[size(image,2) size(image,1), size(image,3)]);
                h5write(fileName,dsetName,permute(image,[2 1 3]));
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
