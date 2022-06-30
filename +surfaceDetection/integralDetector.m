classdef integralDetector < surfaceDetection.surfaceDetector
    % Segmentation based on prediction maps from ilastik that maximize the
    % enclosed probability within a contiguous volume of pixels, given an
    % initial guess. A Chan-Vese energy functional is minimized, with three
    % terms: surface tension, pressure, and attachment energy which fixes
    % the boundary to high gradients in probabilities from ilastik
    % training. Note: bitdepth of output h5 is matched to input image
    % bitdepth.
    % This function uses activecontour() in MATLAB, not morphsnakes.
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
    
    % Copyright 2019-2020 Noah Mitchell
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
        defaultOptions = struct('channel', 1, ...  % which channel to use for integration in Chan-Vese energy functional (if there are only two channels in training, then this does NOT matter)
            'ssfactor', 4,... % subsampling factor: downsampling of raw data
            'niter', 100, ... % how many iterations before exit if no convergence
            'niter0', 100, ... % how many iterations before exit if no convergence for first timepoint
            'pre_pressure', -5, ... % number of dilation/erosion passes for positive/negative values
            'pressure', 0., ... % float: how many pressure (dilation/erosion) steps per iteration
            'tension', 0,... % float: how many smoothing/surface tension steps per iteration (can be <1)
            'post_pressure', 0, ... % how many iterations to dilate (if positive) or erode (if negative) after convergence
            'exit_thres', 1e-6, ... % convergence threshold: maximum difference between subsequent level sets upon which to exit algorithm ('close enough')
            'foreGroundChannel',2, ... % the index of the first dimension of the 4d input data (if 4d)
            'fileName',[], ... % the filename of h5 to train on
            'meshDir', './mesh_output/', ...  % the directory for all output data/images
            'dataDir', './', ...  % the directory for all h5 training data/images
            'ofn_ls', 'ls_%06d.mat', ...  % the output filename for level sets
            'ofn_ply', 'mesh_ls_%06d.ply', ... % the output filename for PLY files
            'timepoint', 0, ... % which timepoint in the data to consider
            'zdim',2, ... % Which dimension is the z dimension
            'ofn_smoothply', 'mesh_%06d.ply',... % the output file name (not including path directory)
            'init_ls_fn', 'none', ... % the name of the initial level set to load, if any
            'run_full_dataset', false, ... % run MS on a time series, not just one file
            'radius_guess', 10, ... % radius of the initial guess sphere
            'dset_name', 'exported_data', ... % the name of the dataset to load from h5
            'center_guess', 'empty_string', ... % xyz of the initial guess sphere ;
            'plot_mesh3d', false, ...  % if save is true, plot intermediate results in 3d 
            'mask', 'none', ... % filename for mask to apply before running MS
            'mesh_from_pointcloud', false, ... % use a pointcloud from the marching cubes algorithm rather than a mesh to create smoothed mesh
            'prob_searchstr', '_Probabilities.h5', ... % if dataset mode, what string to seek for loading all probabilities in data directory (glob datadir/*searchstr)
            'physicalaxisorder', 'xyzc', ... % axis order relative to mesh axis order by which to process the point cloud prediction. To keep as mesh coords, use xyzc
            'preilastikaxisorder', 'xyzc', ... % axis order as output by ilastik probabilities h5. To keep as saved coords use xyzc
            'ilastikaxisorder', 'cxyz', ... % axis order as output by ilastik probabilities h5. To keep as saved coords use xyzc
            'include_boundary_faces', true,... % keep faces along the boundaries of the data volume if true
            'smooth_with_matlab', 0.01, ... % if <0, use meshlab. If >0, smooth the mesh after marching cubes mesh creation using matlab instead of mlxprogram, with diffusion parameter lambda = this value. If =0, no smoothing.
            'target_edgelength', 6, ...             % float, if using Matlab smoothing, target edge length for mesh resampling
            'enforceSingleComponent', false, ...  % enforce that the resulting mesh is a single component
            'enforceQuality', false, ...    % enforce sphere-like topology of output mesh
            'enforceTopology', false, ...       % Enforce that the output mesh must have targetEulerCharacteristic
            'targetEulerCharacteristic', 2, ... % if enforceTopology, what EulerCharacteristic must the output mesh have
            'maxIterRelaxMeshSpikes', 1, ...  % if smooth_with_matlab=true, #iterations to relax cone-like features in the mesh
            'overwrite', false) ;           % overwrite results on disk if they already exist
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
            % none -- simply uses this.options
            
            opts = this.options ;
            
            % unpack opts
            meshDir = opts.meshDir ;
            dataDir = opts.dataDir ;
            pressure = opts.pressure ;
            tension = opts.tension ;
            pre_pressure = opts.pre_pressure ;
            post_pressure = opts.post_pressure ;
            tar_length = opts.target_edgelength ;
            enforceSingleComponent = opts.enforceSingleComponent ;
            enforceQuality = opts.enforceQuality ;
            maxIterRelaxMeshSpikes = opts.maxIterRelaxMeshSpikes ;
            ofn_ls = opts.ofn_ls ;
            ofn_ply = opts.ofn_ply ;
            ofn_smoothply = opts.ofn_smoothply ;
            init_ls_fn = opts.init_ls_fn ;
            tp = opts.timepoint ;
            radius_guess = opts.radius_guess ;
            center_guess = opts.center_guess ;
            niter = opts.niter ; 
            niter0 = opts.niter0 ;
            ssfactor = opts.ssfactor ;
            smooth_with_matlab = opts.smooth_with_matlab ;
            overwrite = opts.overwrite ; 
            enforceTopology = opts.enforceTopology ;
            targetEulerCharacteristic = opts.targetEulerCharacteristic ;
            
            debugMsg(1, ['integralDetector.detectSurface() : channel='...
                num2str(opts.channel), ...
                ', ssfactor=' num2str(opts.ssfactor)...
                ', fileName =' num2str(opts.fileName)...
                ', foreGroundChannel =' num2str(opts.foreGroundChannel)...
                ', zDim =' num2str(opts.zdim),'\n']);
            
            %% load the exported data out of the ilastik prediction
            fn_tmp = sprintf(opts.fileName, tp) ;
            if strcmpi(fn_tmp(end-2:end), '.h5')
                fn_tmp = fn_tmp(1:end-4) ;
            end
            fileName = fullfile(dataDir, ...
                [fn_tmp, '_Probabilities.h5']) ;
            disp(['Searching for init h5 file: ' fileName])
            h5fileInfo = h5info(fileName);
            if strcmp(h5fileInfo.Datasets.Name,'exported_data')
                file = h5read(fileName,'/exported_data');
            elseif strcmp(h5fileInfo.Datasets.Name,'volume')
                file = h5read(fileName,'/volume/prediction');
            else
                error(['Please provide a regular prediction from ilastik, either in', ...
                    'the format of version 1.1 or 0.5 (ie with exported_data as a dataset)']);
            end

            % ilastik internally swaps axes. 1:x, 2:y, 3:z, 4:class
            % strategy: put into xyzc format, then pop last index
            if strcmp(opts.ilastikaxisorder, 'xyzc')
                pred = file ;
            elseif strcmp(opts.ilastikaxisorder, 'yxzc')
                % to convert yxzc to xyzc, put x=2 y=1 z=3 c=4
                pred = permute(file,[2,1,3,4]);
            elseif strcmp(opts.ilastikaxisorder, 'zyxc')
                % to convert yxzc to xyzc, put x=3 y=2 z=1 c=4
                pred = permute(file,[3,2,1,4]);
            elseif strcmp(opts.ilastikaxisorder, 'yzcx')
                % to convert yxzc to xyzc, put x=4 y=1 z=2 c=3
                pred = permute(file,[4,1,2,3]);
            elseif strcmp(opts.ilastikaxisorder, 'cxyz')
                % to convert yxzc to xyzc, put x=2 y=3 z=4 c=1
                pred = permute(file,[2,3,4,1]);
            elseif strcmp(opts.ilastikaxisorder, 'cyxz')
                % to convert yxzc to xyzc, put x=2 y=3 z=4 c=1
                pred = permute(file,[3,2,4,1]);
            elseif strcmp(opts.ilastikaxisorder, 'czyx')
                % to convert yxzc to xyzc, put x=1>4 y=2>3 z=3>2 c=4>1
                pred = permute(file,[4,3,2,1]);
            elseif strcmp(opts.ilastikaxisorder, 'cyzx')
                % to convert cyzx to xyzc put x=1>4 y=2>2 z=3>3 c=4>1
                pred = permute(file,[4,2,3,1]);
            elseif strcmp(opts.ilastikaxisorder, 'cxzy')
                % to convert cxzy to xyzc put x=1>2 y=2>4 z=3>3 c=4>1
                pred = permute(file,[2,4,3,1]);
            else
                error('Have not coded for this axisorder. Do so here')
            end
            
            try
                assert(all(size(pred(:,:,:,1)) > 3))
            catch
                msg = ['One of the spatial dimensions has been labeled as ' ...
                    'a color channel! Change ilastikaxisorder=' ...
                    opts.ilastikaxisorder ' from file size=[' num2str(size(file)) ...
                    '] which translated to [' num2str(size(pred)) ']'] ;
                error(msg) 
            end
            
            %% Define the mesh we seek
            outputLSfn = fullfile(meshDir, sprintf(ofn_ls, tp)) ;
            outputMesh = fullfile(meshDir, sprintf(ofn_ply, tp)) ;

            disp(['does outputMesh exist: ', num2str(exist(outputMesh, 'file'))])
            disp(outputMesh)
           
            if ~exist(outputMesh, 'file')
                if isempty(opts.fileName)
                    error('Please provide a regular prediction from ilastik in h5 format.');
                end

                %---------------------------------
                % Segmentation of a prediction map from ilastik.
                %---------------------------------

                % Check if previous time point's level set exists to use as a seed
                % First look for supplied fn from detectOptions.
                % If not supplied (ie init_ls_fn is none or empty string, then
                % seek previous timepoint output from MS algorithm.

                disp(['init_ls_fn = ', init_ls_fn])
                disp(['ofn_ls = ', ofn_ls])
                if strcmp(init_ls_fn, 'none') || strcmp(init_ls_fn, '')
                    % User has NOT supplied fn from detectOptions
                    init_ls_fn = sprintf(ofn_ls, (tp - 1) ) ;
                end

                % Obtain guess or instructions for where to (1) start the 
                % initial guess level set sphere, and (2) to view the
                % result
                if isempty(center_guess) || strcmpi(center_guess, 'empty_string')
                    centers = size(pred) * 0.5 ;
                    centers = centers(1:3) ;
                else
                    centers = str2num(center_guess) ;
                end

                if ~exist(outputLSfn, 'file') || overwrite
                    disp([ 'initial level set fn = ', init_ls_fn])
                    if exist(init_ls_fn, 'file') || exist([init_ls_fn '.mat'], 'file') 
                        % It does exist, and the given name is the RELATIVE path.    
                        % Use it as a seed (initial level set) 
                        disp('running using initial level set')
                        init_ls = load(init_ls_fn, 'BW') ;
                        init_ls = init_ls.BW ;
                        niter_ii = niter ;
                    elseif exist(fullfile(meshDir, init_ls_fn), 'file') || ...
                            exist(fullfile(meshDir,[init_ls_fn '.mat']), 'file') 
                        % It does exist, and given name is the relative path
                        % without the extension. 
                        % Use it as a seed (initial level set)
                        disp('running using initial level set')
                        init_ls = load(fullfile(meshDir, init_ls_fn), 'BW') ;
                        init_ls = init_ls.BW ;
                        niter_ii = niter ; 
                    else
                        % The guess for the initial levelset does NOT exist, so use
                        % a sphere for the guess.
                        disp(['Using default sphere of radius ' num2str(radius_guess) ...
                            ' for init_ls -- no such file on disk: ' ...
                            fullfile(meshDir, init_ls_fn )])
                        init_ls = zeros(size(squeeze(pred(:, :, :, opts.foreGroundChannel)))) ;
                        SE = strel("sphere", radius_guess) ;
                        SE = SE.Neighborhood ;
                        se0 = size(SE, 1) ;
                        rad = ceil(se0*0.5) ;
                        assert( all(size(SE) == se0)) ;
                        dd = centers - rad ;
                        xmin = max(1, dd(1)) ;
                        ymin = max(1, dd(2)) ;
                        zmin = max(1, dd(3)) ;
                        % xmax = min(size(init_ls, 1), dd(1)+se0);
                        % ymax = min(size(init_ls, 2), dd(2)+se0) ;
                        % zmax = min(size(init_ls, 3), dd(3)+se0) ;

                        % distance from edge of data volume to edge of strel:
                        sz0 = size(init_ls) ;

                        % check if SE is entirely contained within init_ls
                        if all(dd > 0)
                            % minima are within the boundary
                            init_ls(xmin:xmin+se0-1, ymin:ymin+se0-1, ...
                                zmin:zmin+se0-1) = SE ;
                            init_ls = init_ls(1:sz0(1), 1:sz0(2), 1:sz0(3)) ;
                        elseif all(dd + se0 < 0)
                            % the maxima are all contained within the volume
                            init_ls(xmin:xmin+se0+dd(1)-1, ...
                                ymin:ymin+se0+dd(2)-1, ...
                                zmin:zmin+se0+dd(3)-1) = ...
                                SE(-dd(1):se0, -dd(2):se0, -dd(3):se0) ;
                            init_ls = init_ls(1:sz0(1), 1:sz0(2), 1:sz0(3)) ;
                        else
                            % the maxima are all contained within the volume
                            chopx = dd(1) < 0 ;
                            chopy = dd(2) < 0 ;
                            chopz = dd(3) < 0 ;
                            init_ls(xmin:xmin+se0+(dd(1)*chopx)-1, ...
                                ymin:ymin+se0+dd(2)*chopy-1, ...
                                zmin:zmin+se0+dd(3)*chopz-1) = ...
                                    SE(-dd(1)*chopx+1:se0, ...
                                    -dd(2)*chopy+1:se0, ...
                                    -dd(3)*chopz+1:se0) ;
                            init_ls = init_ls(1:sz0(1), 1:sz0(2), 1:sz0(3)) ;
                        end

                        try
                            assert(any(init_ls(:)))
                        catch
                            error('The initial guess is outside the data volume')
                        end

                        niter_ii = niter0 ;
                    end

                    % Flip axis order LR of output mesh to return to MATLAB
                    % orientation? NO, not helpful to do "-permute_mesh 'zyx'"
                    % since already did morphsnakesaxisorder = fliplr() earlier
                    % command = [command ' -adjust_for_MATLAB_indexing'] ;


                    %% Extract contour/isosurface of levelset

                    % Pre-processing
                    if pre_pressure < 0
                        disp(['eroding input LS by pre_pressure=' num2str(pre_pressure)])
                        SE = strel('sphere', abs(pre_pressure)) ;
                        init_ls = imerode(init_ls, SE) ;
                    elseif pre_pressure > 0                
                        disp(['dilating input LS by pre_pressure=' num2str(abs(pre_pressure)) ])
                        SE = strel('sphere', abs(pre_pressure)) ;
                        init_ls = imdilate(init_ls, SE) ;
                    end

                    data = pred(:, :, :, opts.foreGroundChannel) ;
                    % data_clipped = data - 0.1 ;
                    % data_clipped(data_clipped < 0) = 0. ;

                    disp(['niter is ', num2str(niter_ii)]);                
                    BW = activecontour(data, init_ls, niter_ii, 'Chan-Vese', ...
                        'SmoothFactor', tension, 'ContractionBias', -pressure) ;

                    % Post processing
                    if post_pressure < 0
                        disp('eroding result by post_pressure...')
                        SE = strel('sphere', abs(post_pressure)) ;
                        BW = imerode(BW, SE) ;
                    elseif post_pressure > 0          
                        disp('dilating result by post_pressure...')      
                        SE = strel('sphere', abs(post_pressure)) ;
                        BW = imdilate(BW, SE) ;
                    end

                    % preview current results
                    debugLevel = getpref('ImSAnE', 'msgLevel');
                    if debugLevel > 0
                        clf
                        if centers(3) < size(BW, 3) && centers(3) > 0.5 
                            subplot(1, 3, 1)
                            bwPage = squeeze(BW(round(centers(1)), :, :)) ;
                            datPage = squeeze(data(round(centers(1)), :, :)) ;
                            rgb = cat(3, bwPage, datPage, datPage) ;
                            imshow(rgb)
                            subplot(1, 3, 2)
                            bwPage = squeeze(BW(:, round(centers(2)), :)) ;
                            datPage = squeeze(data(:, round(centers(2)), :)) ;
                            rgb = cat(3, bwPage, datPage, datPage) ;
                            imshow(rgb)
                            subplot(1, 3, 3)
                            bwPage = squeeze(BW(:, :,round(centers(3)))) ;
                            datPage = squeeze(data(:,:,round(centers(3)))) ;
                            rgb = cat(3, bwPage, datPage, datPage) ;
                            imshow(rgb)
                            sgtitle('level set found...')

                            pause(1)

                            sgtitle('level set found...')
                            for qq =1:max(size(BW))
                                page = min(qq, size(BW,1)) ;
                                subplot(1, 3, 1)
                                bwPage = squeeze(BW(page, :, :)) ;
                                datPage = squeeze(data(page, :, :)) ;
                                rgb = cat(3, bwPage, datPage, datPage) ;
                                imshow(rgb)
                                page = min(qq, size(BW,2)) ;
                                subplot(1, 3, 2)
                                bwPage = squeeze(BW(:, page, :)) ;
                                datPage = squeeze(data(:, page, :)) ;
                                rgb = cat(3, bwPage, datPage, datPage) ;
                                imshow(rgb)
                                page = min(qq, size(BW,3)) ;
                                subplot(1, 3, 3)
                                bwPage = squeeze(BW(:, :,page)) ;
                                datPage = squeeze(data(:,:,page)) ;
                                rgb = cat(3, bwPage, datPage, datPage) ;
                                imshow(rgb)
                                pause(0.00001)
                            end
                        else
                            xframe = round(size(BW, 1) * 0.5) ;
                            yframe = round(size(BW, 2) * 0.5) ;
                            zframe = round(size(BW, 3) * 0.5) ;
                            subplot(1, 3, 1)
                            bwPage = squeeze(BW(xframe,:,:)) ;
                            datPage = squeeze(data(xframe,:,:)) ;
                            rgb = cat(3, bwPage, datPage, datPage) ;
                            imshow(rgb)
                            subplot(1, 3, 2)
                            bwPage = squeeze(BW(:, yframe, :)) ;
                            datPage = squeeze(data(:, yframe, :)) ;
                            rgb = cat(3, bwPage, datPage, datPage) ;
                            imshow(rgb)
                            subplot(1, 3, 3)
                            bwPage = squeeze(BW(:, :,zframe)) ;
                            datPage = squeeze(data(:,:,zframe)) ;
                            rgb = cat(3, bwPage, datPage, datPage) ;
                            imshow(rgb)
                            sgtitle('level set found...')

                        end
                        % pause to draw the figure and show it in foreground
                        pause(1e-5)
                    end

                    % Remove all but biggest component
                    if enforceSingleComponent
                        CC = bwconncomp(BW, 6);
                        numPixels = cellfun(@numel,CC.PixelIdxList);
                        [~,idx] = max(numPixels);
                        BW = false(size(BW));
                        BW(CC.PixelIdxList{idx}) = true;

                        % bwareaopen(BW, PP, 6) ; % choose connectivity to be six so that face junctions are required
                    end
                else
                    load(outputLSfn, 'BW')
                end

                % Extract mesh from BW
                if opts.include_boundary_faces 
                    % Pad the walls with zeros
                    BW2 = zeros(size(BW) + 2) ;
                    BW2(2:end-1, 2:end-1, 2:end-1) = BW ;

                    % Convert BW to mesh
                    mesh = isosurface(BW2, 0.5) ;
                    mesh.vertices = mesh.vertices - 1 ;
                else
                    mesh = isosurface(BW, 0.5) ;
                end
                
                mesh.vertices = mesh.vertices * ssfactor ;
                % Swap X<->Y axes since MATLAB did this in isosurface/contour
                mesh.vertices = mesh.vertices(:, [2, 1, 3]) ;

                % Check it
                debugLevel = getpref('ImSAnE', 'msgLevel');
                if debugLevel > 2
                    trisurf(triangulation(mesh.faces, mesh.vertices), 'edgecolor', 'none')
                end
                
                % Write init_ls for next timepoint to disk
                save(outputLSfn, 'BW')

                % Write mesh to disk
                disp(['Writing PLY mesh to disk: ' outputMesh])
                plywrite(outputMesh, mesh.faces, mesh.vertices)
            else
                disp(['output PLY already exists: ', outputMesh])
            end

            %% Clean up mesh file for this timepoint using MeshLab/MATLAB 

            % Here use the boundary mesh from marching cubes or MATLAB's 
            % isosurface() to make a smooth mesh
            rawMesh = outputMesh ; 
            outputMesh = fullfile(meshDir, sprintf(ofn_smoothply, tp));

            disp(['outputMesh = ', outputMesh])
            %bad = so_bad
            if ~exist( outputMesh, 'file') || overwrite
                % Smooth with either meshlab or matlab
                if smooth_with_matlab < 0
                    % USE MESHLAB, not matlab
                    if use_pointcloud
                        % Use the pointcloud from the level set rather than the
                        % boundary mesh from marching cubes
                        %----------------------------------------------------------------------
                        % Extract the implicit level set as a 3D binary array
                        %----------------------------------------------------------------------

                        % The file name of the current time point
                        % The 3D binay array
                        bwLS = h5read( outfnLS, '/implicit_levelset' );

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
                        command = ['meshlabserver -i ' pointCloudFileName, ...
                            ' -o ' outputMesh, ' -s ' mlxprogram ' -om vn'] ;
                        disp(['running ' command])
                        system( command );
                    else
                        % Use the marching cubes mesh surface to smooth
                        command = ['meshlabserver -i ' infile ' -o ' outputMesh, ...
                            ' -s ' mlxprogram ' -om vn'];
                        % Either copy the command to the clipboard
                        clipboard('copy', command);
                        % or else run it on the system
                        disp(['running ' command])
                        system(command)
                    end
                elseif smooth_with_matlab == 0
                    disp('No smoothing, with either matlab or meshlab')
                    mesh = read_ply_mod(infile) ;
                    disp('Compute normals...')
                    mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
                    plywrite_with_normals(outputMesh, mesh.f, mesh.v, mesh.vn)
                elseif smooth_with_matlab > 0 
                    % Smooth with MATLAB
                    disp(['Smoothing with MATLAB using lambda = given value of ' num2str(smooth_with_matlab)])
                    mesh = read_ply_mod(rawMesh) ;

                    % Check that this behaves the way we want
                    % mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
                    % size(mesh.vn)
                    % size(mesh.v)
                    % assert(size(mesh.vn, 1) == size(mesh.v, 1))

                    if enforceSingleComponent
                        [ F, V, oldVertexIDx, C ] = ...
                            remove_isolated_mesh_components(mesh.f, mesh.v) ;
                    else
                        F = mesh.f ;
                        V = mesh.v ;
                    end

                    num_iter = 5;                
                    protect_constraints = false;
                    sm_opts = struct('lambda', smooth_with_matlab, ...
                        'tar_length', tar_length, ...
                        'num_iter', num_iter, ...
                        'protect_constraints', protect_constraints, ...
                        'enforceQuality', enforceQuality, ...
                        'maxIterRelaxMeshSpikes', maxIterRelaxMeshSpikes, ...
                        'enforceTopology', enforceTopology, ... 
                        'targetEulerCharacteristic', targetEulerCharacteristic ) ;
                    [V, F] = remesh_smooth_iterate(V,F, sm_opts) ;

                    % newV = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], smooth_with_matlab) ;
                    disp('Compute normals...')
                    mesh.v = V ;
                    mesh.f = F ;
                    mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
                    disp(['Saving smoothed mesh to ' outputMesh])
                    plywrite_with_normals(outputMesh, mesh.f, mesh.v, mesh.vn)
                end
            else
                disp(['t=', num2str(tp) ': smoothed mesh file found on disk...'])
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if exist(outputMesh, 'file')
                disp(['reading PLY ', outputMesh])
                tmp = read_ply_mod(outputMesh);
                vv = tmp.v ;
                points = struct('x', vv(:, 1), 'y', vv(:, 2), 'z', vv(:, 3));

                % Unlike in other detectors, here we have already rescaled
                % the data
                x = cat(1,points.x);
                y = cat(1,points.y);
                z = cat(1,points.z);
                pointCloud = [x,y,z];

                %--------------------------------------------------
                % scale point cloud to full size and set alignment
                %--------------------------------------------------
                largePC = zeros(size(pointCloud));
                for i = 1:3
                    largePC(:,i) = (pointCloud(:,i)-1) + 1;
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
                xSize = size(pred, 1);
                ySize = size(pred, 2);
                zSize = size(pred, 3);
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
            % Prepare stack for Ilastik segmetnation. This outputs an h5
            % file of subsampled intensity data on which to train.
            %
            % prepareIlastik(stack)
            
            % Accoring to the specified options sub-sample the stack and
            % save it for analysis with ilastik.
            
            opts = this.options;
            
            im = stack.image.apply();
            
            if ~strcmpi(opts.fileName(end-2:end), '.h5')
                fileName = [opts.fileName,'.h5'];
            else
                fileName = opts.fileName ;
            end
            
            if exist(fileName,'file')
                delete(fileName)
            end
            
            dsetName = '/inputData';
            
            % Determine the axis order for ilastik training data
            if strcmp(opts.preilastikaxisorder, 'xyzc') 
                axperm = [1 2 3 4] ;
            elseif strcmp(opts.preilastikaxisorder, 'yxzc')
                axperm = [2 1 3 4] ;
            elseif strcmp(opts.preilastikaxisorder, 'zxyc')
                axperm = [3 1 2 4] ;
            elseif strcmp(opts.preilastikaxisorder, 'czxy')
                axperm = [4 3 1 2] ;
            elseif strcmp(opts.preilastikaxisorder, 'czyx')
                axperm = [4 3 2 1] ;
            elseif strcmp(opts.preilastikaxisorder, 'cxyz')
                axperm = [4 1 2 3] ;
            else
                error(['Have not coded for this axis permutation yet: ', ...
                    opts.preilastikaxisorder])
            end
            
            % Subsample the image to save for ilastik training
            for c = 1 : length(im)
                image(:,:,:,c) = im{c}(1:opts.ssfactor:end,1:opts.ssfactor:end,1:opts.ssfactor:end);
            end
            
            % Now save the subsampled images to h5 using the axis order
            % specified by axperm
            if ndims(image)==4
                disp(['Writing file: ' fileName])
                if isa(image, 'uint8')
                    h5create(fileName,dsetName,[size(image,axperm(1)) size(image,axperm(2)), ...
                        size(image,axperm(3)) size(image,axperm(4))], ...
                        'datatype', 'uint8',...
                        'Chunksize', [size(image,axperm(1)) size(image,axperm(2)),...
                        size(image,axperm(3)) size(image,axperm(4))]);
                elseif isa(image, 'uint16')
                    h5create(fileName,dsetName,[size(image,axperm(1)) size(image,axperm(2)),...
                        size(image,axperm(3)) size(image,axperm(4))], ...
                        'Datatype', 'uint16');
                else
                    error('Did not recognize bitdepth of image. Add capability here')
                end
                h5write(fileName,dsetName,permute(image, axperm));
            else
                % truncate the axis permutation to include just 3 dims
                axperm = axperm(1:3) ;
                disp(['Writing file: ' fileName])
                
                if isa(image, 'uint8')
                    h5create(fileName,dsetName,[size(image,axperm(1)) size(image,axperm(2)),...
                        size(image,axperm(3))], ...
                        'datatype', 'uint8',...
                        'Chunksize', [size(image,axperm(1)) size(image,axperm(2)), size(image,axperm(3)) ]) ; 
                elseif isa(image, 'uint16')
                    h5create(fileName,dsetName,...
                        [size(image,axperm(1)) size(image,axperm(2)), size(image,axperm(3))], ...
                        'Datatype', 'uint16')
                else
                    error('Did not recognize bitdepth of image. Add capability here')
                end
                h5write(fileName,dsetName,permute(image, axperm));
            end
            
        end
        
        % ------------------------------------------------------
        % Check point cloud against original image
        % ------------------------------------------------------
        function inspectQuality(this, inspectOpts, stack)
            %   inspect quality of fit in single slice in dimension specified
            %   by options and display image.
            %   NOTE: Unlike in other detectors, here the mesh is already
            %   scaled by ssfactor, so no need to override inspectQuality 
            %   to deal with subsampled point cloud.
            %
            %   inspectQuality(inspectOpts, stack)
            %
                        
            inspectQuality@surfaceDetection.surfaceDetector(this, inspectOpts, stack);
        end
        
    end
end
