function [geodesicPaths, pointLocations] = surfaceGeodesicPairs( face, ...
    vertex, pointPairs, pointCoordinates )
%SURFACEGEODESICPAIRS This function calculates the surface geodesic paths
%between pairs of points on a 3D mesh triangulation.  The points are
%assumed to be elements of face interiors, i.e. not vertices or elements of
%triangulation edges.  Geodesic paths are returned as a series of 3D
%points.
%   Detailed explanation goes here

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------

if(nargin<1), error('Please supply face connectivity list!'); end
if(nargin<2), error('Please supply vertex coordinates!'); end
if(nargin<3), error('Please supply source/target coordinates!' ); end
if(nargin<4), error('Please supply source/target pair list!' ); end

% Check the size of the face connectivity list
sizef = size(face);
if ((sizef(2)~=3)||(length(sizef)~=2))
    error('The face list is improperly sized');
end

% Check the size of the vertex list
sizev = size(vertex);
if ((sizev(2)~=3)||(length(sizev)~=2))
    error('The vertex list is improperly sized');
end

% Check if vertex indices exist
if ( max(face(:)) > sizev(1) )
    error('The face list contains an undefined vertex');
elseif ( min(face(:)) < 1 )
    error('The face list contains a vertex index smaller than 1');
end

% Check the size of the point location list
sizepc = size(pointCoordinates);
if ((sizepc(2)~=3)||(length(sizepc)~=2))
    error('The source/target coordinate list is improperly sized');
end

% Check the size of the pair list
sizepp = size(pointPairs);
if ((sizepp(2)~=2)||(length(sizepp)~=2))
    error('The source/target pair list is improperly sized');
end

% Check if point indices exist
if ( max(pointPairs(:)) > sizepc(1) )
    error('The source/target pair list contains an undefined point');
elseif ( min(pointPairs(:)) < 1 )
    error('The source/target pair list contains a point smaller than 1');
end

% *************************************************************************
% TODO: CREATE BOND LIST SORT TO MINIMIZE SEQUENCE TREE RECONSTRUCTIONS IN
% NUMERICAL ALGORITHM
% *************************************************************************

% Update vertex and point IDs to match the 0-indexing in C++
face = face-1; pointPairs = pointPairs-1;

%--------------------------------------------------------------------------
% Calculate Geodesic Paths!
%--------------------------------------------------------------------------

[ geodesicPaths, ...
    pointFaceLocations, ...
    pointBarycentricCoordinates ] = surface_geodesic_pairs( uint64(face), ...
    vertex, uint64(pointPairs), pointCoordinates );

%--------------------------------------------------------------------------
% Output Processing
%--------------------------------------------------------------------------

% Add +1 to accound for MATLAB 1-indexing
pointFaceLocations = pointFaceLocations + 1;

% *************************************************************************
% TODO: CREATE MAP FROM SORTED BOND LIST TO INPUT BOND LIST
% *************************************************************************

% Assemble cell location struct -------------------------------------------

pfl_Cell = cell( sizepc(1), 1 );
pbc_Cell = cell( sizepc(1), 1 );
for i = 1:sizepc(1)
    
    pfl_Cell{i} = pointFaceLocations(i);
    pbc_Cell{i} = pointBarycentricCoordinates(i,:);
    
end

pointLocations = struct( 'face', pfl_Cell, ...
    'barycentricCoordinates', pbc_Cell );

end


