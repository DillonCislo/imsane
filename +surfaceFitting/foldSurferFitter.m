classdef foldSurferFitter < surfaceFitting.surfaceFitter
    % Generate a SurfaceOfInterest object with atlas containing  
    % charts from an externally produced triangular mesh representation of 
    % a surface with initially spherical topology. 
    % This is suitable for singly-connected meshes that can be
    % parameterized by a cetnerline curve and a folded (multivalued) net of
    % 2D rectilinear coordinates.
       
    %---------------------------------------------------------------------
    % license
    %---------------------------------------------------------------------

    % Copyright 2022 DJCislo and NPMitchell
    % extension to: ImSAnE by Idse Heemskerk and Sebastian Streichan
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
    properties (SetAccess = protected)
        
        % charts - charts this fitter can produce and their properties
        %
        % the rows have the structure:
        % {name, description, stepsize, intersections, fundamental, desired}
        %
        % stepsize:         2d vector of stepsize
        % intersections:    lists indices of charts with overlapping domain
        % fundamental:      boolean, does chart define a set in the topology
        % desired:          boolean indicating whether to produce it
        charts = struct( ...
            'name', {'uv'},...
            'description', {'Dirichlet mapping to the unit square'},...
            'stepSize', {[1 1]},...
            'intersections', {[]},...
            'fundamental', {1},...
            'desired', {1});
        
        
        foldSurfer 
        
        
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = foldSurferFitter(xpMeta, foldSurferOpts)
            % MESHWRAPPER Create a SOI from a TubULAR compatible mesh
            %
            % meshWrapper(xpMeta, opts)
            %
            % xpMeta: ImSAnE experiment meta data
            % opts: TubULAR instantiation options
            %

            if (nargin < 2)
                foldSurferOpts = struct() ;
                % error(['You must pass the metadata of the current ' ...
                %     'experiment to foldSurferFitter as: ' ...
                %     'xp = struct(fileMeta, expMeta, detectOptions)']);
            end
                                    
            % call superclass constructor            
            this = this@surfaceFitting.surfaceFitter();
                        
            % initialize fitOptions
            this.fitOptions.fixAxis         = 1;
            this.fitOptions.rotation        = eye(4);
            this.fitOptions.translation     = eye(4);
            this.fitOptions.fixResolution   = 0;
            this.fitOptions.resolution      = [];
            this.fitOptions.phase           = 0;
            
            % initialize fittedParam
            this.fittedParam = struct() ;
            this.fittedParam.mesh = struct() ;
            this.fittedParam.im = struct() ;
            
        end
        
        % ------------------------------------------------------
        % fitting
        % ------------------------------------------------------
        
        function fitSurface(this, mesh, options)
            % FITSURFACE Convert meshes generated using TubULAR into an
            % ImSAnE-style surface representation
            %
            % fitSurface(mesh)
            %
            % mesh: struct with fields:
            %       v:    Nx3 vertex positions
            %       vn:   Nx3 vertex normals
            %       f:    Mx3 faces (triangles)
            % options : optional struct with fields
            %       options : 
            %
            % fitOptions can be set through the generic function
            % setFitOptions. The tubularFitter defines the following:
            %
            % see also surfaceFitting.surfaceFitter.setFitOptions
            preview = false ;
            if nargin < 3
                options = struct() ;
            else
                if isfield(options, 'preview')
                    preview = options.preview ;
                end
            end
            
            % Add current timepoint's mesh to fittedParam            
            v3d = mesh.v;
            FF = mesh.f;
            vn = mesh.vn;
            
            % make sure normals are normalized
            normalnorm = sqrt(vn(:,1).^2 + vn(:,2).^2 + vn(:,3).^2);
            for i=1:3
                vn(:,i) = vn(:,i)./normalnorm;
            end
            mesh.vn = vn;
            
            % store the mesh in fittedParam
            this.fittedParam.mesh = mesh;
            
            % fittedPoints are mesh vertices for now
            this.fittedPoints = {mesh.v(:,1), mesh.v(:,2), mesh.v(:,3)};
            
            %--------------------------------------------------------------
            % generate submeshes for patches on which charts will be defined
            %--------------------------------------------------------------   
            if isfield(mesh, 'cornerIDx')
                cornerIDx = mesh.cornerIDx ;
            else
                % try to find the corners
                TR3D = triangulation(FF, v3d);
                E = TR3D.edges; % #Ex2 edge connectivity list
                bdyIDx = TR3D.freeBoundary; % #BEx2 sorted list of boundary edges

                minx = min(mesh.v(bdyIDx, 1)) ;
                miny = min(mesh.v(bdyIDx, 2)) ;
                maxx = max(mesh.v(bdyIDx, 1)) ;
                maxy = max(mesh.v(bdyIDx, 2)) ;
                rr = vecnorm(mesh.v(bdyIDx, [1,2]), 2, 2) ;
                Ly = maxy - miny ;
                Lx = maxx - minx ;
                candx = bdyIDx(v3d(bdyIDx, 2) < miny+Ly*0.01) ;
                candy = bdyIDx(v3d(bdyIDx, 1) < minx+Lx*0.01) ;
                [~,idx] = max(v3d(candx, 1)) ;
                [~,idy] = max(v3d(candy, 2)) ;
                idx = candx(idx) ;
                idy = candy(idy) ;
                [~,id0] = min(rr) ;
                [~,idr] = max(rr) ;
                id0 = bdyIDx(id0) ;
                idr = bdyIDx(idr) ;
                cornerIDx = [idy, id0, idx, idr] ; % if blank, use default corners
            
                % Check result
                if preview
                    trisurf(triangulation(FF, v3d), 'edgecolor', 'none') ;
                    hold on;
                    for ii = 1:4
                        scatter3(v3d(cornerIDx(ii), 1), v3d(cornerIDx(ii), 2), ...
                            v3d(cornerIDx(ii), 3), 100, 'filled')
                    end
                    axis equal
                    pause(1) ;
                end
            end
            [mesh.u, mesh.cornerIDx] = flattenSquare( mesh.f, mesh.v, cornerIDx ) ;
            
            % Update mesh in fittedParam to include u and cornerIDx
            this.fittedParam.mesh = mesh ;
            
            % Output texturepatch image
            Options.numLayers = [ 15,0 ] ;
            Options.layerSpacing = 0.5 ;
            im = texture_patch_to_image(FF, mesh.u, FF, ...
                v3d(:, [2, 1, 3]), options.IV, Options) ;
            tmp = max(im(:, :, :, 8:12), [], 4) ;
            imshow(tmp) ;
            
            this.fittedParam.im = tmp ;
        end
        
        
        % ------------------------------------------------------
        %  generate embedding
        % ------------------------------------------------------
               
        function [embedding, chart] = generateEmbedding(this, chartName)
            % GENERATEEMBEDDING Create the embedding in some chart
            %
            % [embedding, chart] = generateEmbedding(chartName)
            %
            % The embedding is generated using the result of a fit, stored
            % in fittedParam, so fitSurface should be called first.
            % The available charts are listed in the charts property of the
            % SurfaceFitter object.
            % 
            % See also surfaceFitting.tpsFitter.fitSurface
                        
            % --------------------------------------------
            % Define Embeddings for this chart type
            % --------------------------------------------
            
            debugMsg(2, ['generateEmbedding(' chartName ')\n']);

            % convert to structure used by matlabmesh toolbox
            %subm_tb = makeMesh(subm.v, subm.f, subm.vn);
            %u = embedSCP(subm_tb, 'robust'); %'fiedler', 'generalized'
            % 
            % mesh = this.fittedParam.mesh{ ...
            %     strcmp({this.charts(:).name}, chartName) }

            mesh = this.fittedParam.mesh ;
            
            if strcmpi(chartName, 'uv')
                uv = mesh.u  ;
            else
                error('code for this chart here')
            end
            
            % Interpolate on grid
            u = uv(:,1);
            v = uv(:,2);

            notInf = ~isinf(u) & ~isinf(v);
            u = u(notInf);
            v = v(notInf);

            minu = round(min(u(:)));
            maxu = round(max(u(:)));
            minv = round(min(v(:)));
            maxv = round(max(v(:)));

            uarr = linspace(minu,maxu,500) ;
            varr = linspace(minv,maxv,500) ;
            [ugrid, vgrid] = meshgrid(uarr, varr) ;
            du = uarr(2) - uarr(1) ;
            dv = varr(2) - varr(1) ;
            
            uidx = 1;
            mask = mesh_mask(mesh, uidx, minv, minu, maxv-minv+1, maxu-minu+1);

            grids = {};
            for embIdx = 1:3
                X = mesh.v(notInf, embIdx);
                F = scatteredInterpolant(u,v,X,'linear','none');
                grid = F(ugrid, vgrid);
                grid(~mask) = NaN;
                grids{embIdx} = grid;
            end

            % define the image and domain of the chart
            boundary = {[minu, maxu], [minv, maxv]};
            stepSize = [du dv];
            chartim  = diffgeometry.FiniteSet(chartName,...
                                            boundary, stepSize);
            domain   = chartim.makeIndexSet();

            % chart: name_index -> name
            % the domain is the image of the chart, the image is the domain of
            % the fit, i.e. the embedding space
            chart = diffgeometry.CoordinateMap(domain, chartim,...
                                            chartim.makeHandles());

            % for charts made with FiniteSet.makeHandles we can also set an
            % analytic inverse (for speed and accuracy)
            chartInv = diffgeometry.CoordinateMap(chartim, domain,...
                                    chartim.makeInverseHandles());
            chart.setInverse(chartInv);

            % embedding: cylinder1 -> targetSpace
            embedding = diffgeometry.CoordinateMap(chartim,...
                                     this.fitDomain, grids);

            % embedding \circ chart: cylinder1_index -> targetSpace
            embedding = embedding.compose(chart); 

        end
        
        function normallyEvolve(this, shift)
            % NORMALLYEVOLVE Normally shift surface out or in by some amount
            %
            % normallyEvolve(shift)
            %
            % shift : pixel distance along the normal
            mesh = this.fittedParam.mesh ;
            v = mesh.v;
            vn = mesh.vn;
            this.fittedParam.mesh.v = v + shift*vn;
            
        end

        %------------------------------------------------------
        % populate SOI
        %------------------------------------------------------
        
        function  populateSOI(this, SOI, varargin)
            % POPULATESOI Add fit result to SurfaceOfInterest object
            %
            % populateSOI(SOI)
            % populateSOI(SOI, t)
            %
            % SOI:  SurfaceOfInterest object
            % t:    time, needs to be provided if SOI.dynamic = true
            %
            % Generates chart domain, chart and embedding using the result 
            % of a fit, adds these to the SOI.
            % Fit results are stored in fittedParam, so fitSurface should 
            % be called first.
            % 
            % See also surfaceFitting.tpsFitter.fitSurface
            
            % gti : geometric time index (one for static geometry, time index for
            % dynamic geometry)
            if SOI.dynamic == false
                gti = 1;
                if length(varargin) == 1
                    debugMsg(1, 'static SOI: ignoring time provided\n');
                end
            else
                if length(varargin) == 1
                    gti = SOI.tIdx(varargin{1});
                else
                    error('for dynamic SOI, time argument needs to be provided');
                end
            end
            
            %-----------------------------------------------------------------
            % Define uv / sphi / sphi_sm / ricci 
            %-----------------------------------------------------------------

            for ii = 1:length(this.charts)
                name = this.charts(ii).name ;
                if this.charts(ii).desired
                    
                    debugMsg(2, ['generating ' name ' charts \n']);
                    [embedding, chart] = this.generateEmbedding(name);

                    intersects = {};
                    SOI.topologicalSpace(gti).addSet(chart.domain, intersects);
                    SOI.atlas(gti).addChart(chart);
                    SOI.embedding(gti).addPatch(embedding);
                end
            end
            
           
        end
        
    end
end