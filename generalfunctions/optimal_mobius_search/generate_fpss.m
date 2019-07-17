function [ M, B, eps ] = generate_fpss( movF, movV3, movV2, ...
    tarF, tarV3, tarV2, varargin )
%GENERATE_FPS Generates the 'Farthest Point Sampling Set' for
%Mobius alignment of unit disk domains
%   Input Arguments:
%       - movF:      The connectivity list defining the faces of the moving
%                    mesh
%       - movV3:     The (X,Y,Z) coordinates of the points the 3D mesh
%       - movV2:     The (U,V) coordinates of the points in the 2D mesh
%
%       - tarF:      The connectivity list defining the faces of the target
%                    mesh
%       - tarV3:     The (X,Y,Z) coordinates of the points in the 3D mesh
%       - tarV2:     The (U,V) coordinates of the points in the 2D mesh
%
%   Output Arguments
%       - M:         The set of correspondence points in the moving mesh
%                    Stored as an (Nx1) complex vector
%
%       - B:         The set of correspondence points in the target mesh
%                    Stored as an {NxNj} cell array
%
%       - eps:       The radii of the target circles
%                    Stored as an (Nx1) vector

%==========================================================================
%                           Process Input Data
%==========================================================================

mov3D_b = triangulation( movF, movV3 );
mov2D_b = triangulation( movF, movV2 );

tar3D = triangulation( tarF, tarV3 );
tar2D = triangulation( tarF, tarV2 );

%--------------------------------------------------------------------------
% Construct reduced moving triangulation with the true boundary faces
% removed.  This step is performed so that none of the sample points lie on
% the boundary of the unit disk
%--------------------------------------------------------------------------

movBdy = freeBoundary(mov3D_b); movBdy = unique(movBdy(:));

mov3D = remove_vertex( movBdy, mov3D_b );
mov2D = remove_vertex( movBdy, mov2D_b );

numPnts = 100; % The number of sample points
seedIDx = mov2D.nearestNeighbor([0,0]); % The initial seed point for sampling
method  = 'FPS';

for i = 1:length(varargin)
    if isa(varargin{i}, 'double')
        continue;
    end
    if ~isempty(regexp(varargin{i}, '^[Nn]um[Pp]nts', 'match'))
        numPnts = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Ss]eed[Pp]nt', 'match'))
        seedIDx = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Mm]ethod', 'match'))
        method = varargin{i+1};
    end
    
end

assert(numPnts > 2, 'Correspondence should have at least three points');
assert(numPnts <= size(mov2D.Points,1), ...
    ['Maximum number of correspondence points is: ', ...
    num2str(size(mov2D.Points,1))]);

if size(seedIDx,2) == 2
    seedIDx = mov2D.nearestNeighbor(seedIDx);
elseif size(seedIDx, 2) == 3
    seedIDx = mov3D.nearestNeighbor(seedIDx);
end

assert( size(seedIDx,2) == 1, ...
    ['Seed points must be input as a set of vIDs,' ...
    'a set of (u,v) coordinates, or a set of (x,y,z) coordinates']);

%==========================================================================
%   Determine the remaining sample points
%==========================================================================

if strcmp(method, 'FPS')
    % Farthest Point Sampling 
    
    % vID of each sample point
    sampleIDx = zeros(numPnts,1);
    sampleIDx(1:length(seedIDx)) = seedIDx;
    
    Dmin = inf(size(mov3D.Points,1),length(seedIDx));
    
    % Find the distances for the initial seed points
    for i = 1:length(seedIDx)
        Dmin(:,i) = perform_fast_marching_mesh(mov3D.Points, ...
            mov3D.ConnectivityList,  seedIDx(i));
    end
    
    % The minimum distance of each vertex in the mesh to
    % any of the extant sample points
    Dmin = min(Dmin,[],2);
    
    while i < numPnts
        
        i = i+1;
        
        [~, newID] = max(Dmin);
        assert(~ismember(newID, sampleIDx), ...
            'Update vertex is already a sample point!');
        sampleIDx(i) = newID;
        
        % Update the minimum distances
        D = perform_fast_marching_mesh(mov3D.Points, ...
            mov3D.ConnectivityList,  newID);
        Dmin = min(Dmin, D);
        
    end
    
else
    % Random Point Sampling
    
    sampleIDx = datasample( 1:size(mov2D.Points,1), numPnts, ...
        'Replace', false )';
end

%==========================================================================
%       Identify point sets for Optimal Mobius Alignment
%==========================================================================

% Output Data
eps = zeros(size(sampleIDx));
B   = cell(size(sampleIDx));

M = complex( mov2D.Points(sampleIDx,1), mov2D.Points(sampleIDx,2) );


nn_vIDx = tar3D.nearestNeighbor(mov3D.Points(sampleIDx,:));
centers = tar2D.Points(nn_vIDx,:);
centers = complex(centers(:,1), centers(:,2));

edges = tar2D.edges;

for i = 1:numPnts
    
    % The edges containing the current sample point 
    eIDx = find( any( edges == nn_vIDx(i) , 2 ) );
    uvs = tar2D.Points(edges(eIDx,1),:);
    uve = tar2D.Points(edges(eIDx,2),:);
    
    eps(i) = max(sqrt(sum((uvs - uve).^2,2)));
    B{i} = centers(i);
    
end
    

end

