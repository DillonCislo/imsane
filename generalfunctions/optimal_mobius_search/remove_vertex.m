function [ trNew ] = remove_vertex( vIDx, tr )
%REMOVE_VERTEX Removes specific vertices and any faces containing those
%vertices from a triangulation.
%   Input Parameters:
%       - vIDx:             The list of vertex IDs to be removed
%       - tr:               A triangulation object
%
%   Output Parameters:
%       -trNew:             The updated triangluation

% Remove the vertices at the specified index values
newVertices = tr.Points;
newVertices(vIDx,:) = [];

% Find the new index for each of the remaining vertices
[~, newVertexIndex] = ismember(tr.Points, newVertices, 'rows');

% Find any faces that used the old vertices and remove them
newFaceList = tr.ConnectivityList;
row_v = false(size(newFaceList,1),1);

for i = 1:length(vIDx)
    row_v = row_v | any( newFaceList == vIDx(i), 2 );
end

newFaceList(row_v, :) = [];

% Now update the vertex indices to the new ones
newFaceList = newVertexIndex(newFaceList);

% Remove any unused vertices in the new triangulation
nullVIDx = find( ~ismember( 1:size(newVertices,1), newFaceList ) );
if ~isempty(nullVIDx)
    
    tempVtx = newVertices;
    newVertices(nullVIDx,:) = [];
    [~, newVertexIndex] = ismember(tempVtx, newVertices, 'rows');
    
    newFaceList = newVertexIndex(newFaceList);
    
end

trNew = triangulation(newFaceList, newVertices);


end

