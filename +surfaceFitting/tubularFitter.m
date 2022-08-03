classdef tubularFitter < surfaceFitting.surfaceFitter
    % Generate a SurfaceOfInterest object with atlas containing  
    % charts from an externally produced triangular mesh representation of 
    % a surface with initially spherical topology. 
    % This is suitable for singly-connected meshes that can be
    % parameterized by a cetnerline curve and a tube-like net of
    % longitudinal and circumferential coordinates.
       
    %---------------------------------------------------------------------
    % license
    %---------------------------------------------------------------------

    % Copyright 2022 DJCislo and NPMitchell
    % part of ImSAnE by Idse Heemskerk and Sebastian Streichan
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
            'name', {'uv', 'sphi', 'sphi_sm', 'ricci'},...
            'description', {'Orbifold-Tutte mapping to the unit square, periodic in Y',...
            'The s axis is (circumferentially-averaged) proper length traversed along the longitudinal axis, phi axis is same as v at t=t0, otherwise rotated to match previous timepoint in 3d',...
            'sphi coordinates smoothed over time',...
            'result of a Ricci flow to a rectangle that is periodic in Y'},...
            'stepSize', {[1 .1], [1 .1], [1 1], [1 1]},...
            'intersections', {[],[],[],[]},...
            'fundamental', {1, 1, 0, 0},...
            'desired', {1, 1, 0, 0});
        
        
        tubi 
        
        
    end
    
    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        % ------------------------------------------------------
        % constructor
        % ------------------------------------------------------
        
        function this = tubularFitter(xpMeta, tubularOpts)
            % MESHWRAPPER Create a SOI from a TubULAR compatible mesh
            %
            % meshWrapper(xpMeta, opts)
            %
            % xpMeta: ImSAnE experiment meta data
            % opts: TubULAR instantiation options
            %

            if (nargin < 2)
                error(['You must pass the metadata of the current ' ...
                    'experiment to tubularMeshWrappper as: ' ...
                    'xp = struct(fileMeta, expMeta, detectOptions)']);
            end
                                    
            % call superclass constructor
            this = this@surfaceFitting.surfaceFitter();
            
            % initialize TubULAR object
            this.tubi = TubULAR(xpMeta, tubularOpts) ;
            
            % initialize fitOptions
            this.fitOptions.fixAxis         = 1;
            this.fitOptions.rotation        = eye(4);
            this.fitOptions.translation     = eye(4);
            this.fitOptions.fixResolution   = 0;
            this.fitOptions.resolution      = [];
            this.fitOptions.phase           = 0;
            
            % initialize fittedParam
            this.fittedParam = struct() ;
            % this.fittedParam.meshes = {[],[],[],[]} ;
            % this.fittedParam.cutMeshes = {[],[],[],[]} ;
            this.fittedParam.meshes = {} ;
            this.fittedParam.meshes{1} = [] ;
            this.fittedParam.meshes{2} = [] ;
            this.fittedParam.meshes{3} = [] ;
            this.fittedParam.meshes{4} = [] ;
            this.fittedParam.cutMeshes = {} ;
            this.fittedParam.cutMeshes{1} = [] ;
            this.fittedParam.cutMeshes{2} = [] ;
            this.fittedParam.cutMeshes{3} = [] ;
            this.fittedParam.cutMeshes{4} = [] ;
            
        end
        
        % ------------------------------------------------------
        % fitting
        % ------------------------------------------------------
        
        function fitSurface(this)
            % FITSURFACE Convert meshes generated using TubULAR into an
            % ImSAnE-style surface representation
            %
            % fitSurface()
            %
            % fitOptions can be set through the generic function
            % setFitOptions. The tubularFitter defines the following:
            %
            % see also surfaceFitting.surfaceFitter.setFitOptions
            
            
            % Add current timepoint's mesh to fittedParam
            if this.charts(1).desired
                % uvmesh = this.tubi.getCurrentUVCutMesh() ;
                % uvglue = glueCylinderCutMeshSeam(uvmesh) ;
                % Note here we load the high-resolution (non-downsampled)
                % version of the uv cut mesh. It does not have a
                % rectilinear grid structure to the mesh topology.
                uvmesh = this.tubi.getCurrentCutMesh() ;
                uvglue = glueCylinderCutMeshSeam(uvmesh) ;
            else
                uvmesh = [] ;
                uvglue = [] ;
            end
            if this.charts(2).desired
                sphimesh = this.tubi.getCurrentSPCutMesh() ;
                sphimesh = struct('f', sphimesh.f, ...
                    'v', sphimesh.v, ...
                    'u', sphimesh.sphi, ...
                    'pathPairs', sphimesh.pathPairs) ;
                sphiglue = glueCylinderCutMeshSeam(sphimesh) ;
            else
                sphimesh = [] ;
                sphiglue = [] ;
            end
            if this.charts(3).desired
                sphism_mesh = this.tubi.getCurrentSPCutMeshSm() ;
                sphism_mesh = struct('f', sphism_mesh.f, ...
                    'v', sphism_mesh.v, ...
                    'u', sphism_mesh.u, ...
                    'pathPairs', sphism_mesh.pathPairs) ;
                sphism_glue = glueCylinderCutMeshSeam(sphism_mesh) ;
            else
                sphism_mesh = [] ;
                sphism_glue = [] ;
            end
            if this.charts(4).desired
                tmp = this.tubi.getCurrentRicciMesh() ;
                tmp = tmp.rectangle ;
                ricci_mesh = struct('f', tmp.f, ...
                    'v', tmp.v , ...
                    'u', tmp.u, ...
                    'pathPairs', tmp.pathPairs) ;
                ricci_glue = glueCylinderCutMeshSeam(ricci_mesh) ;
            else
                ricci_mesh = [] ;
                ricci_glue = [] ;
            end
            this.fittedParam.cutMeshes = {uvmesh, sphimesh, sphism_mesh, ricci_mesh} ; 
            this.fittedParam.meshes = {uvglue, sphiglue, sphism_glue, ricci_glue} ;
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
            
            mesh = this.fittedParam.cutMeshes{ ...
                strcmp({this.charts(:).name}, chartName) } ;
            uv = mesh.u  ;

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
        
        function normallyEvolve(this, shift, chartName)
            % NORMALLYEVOLVE Normally shift surface out or in by some amount
            %
            % normallyEvolve(shift)
            %
            % shift : pixel distance along the normal
                        
            chartIdx = strcmp(this.chart.name, chartName) ;
            mesh = this.fittedParam.cutMeshes{chartIdx} ;
            v = mesh.v;
            vn = mesh.vn;
            this.fittedParam.cutMeshes{chartIdx}.v = v + shift*vn;
            
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