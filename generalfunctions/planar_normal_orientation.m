function [ pm ] = planar_normal_orientation( mesh, v2D )
%PLANAR_NORMAL_ORIENTATION Determines whether the prescribed normal vector
%on a mesh triangulation of a surface embedded in 3D corresponds to the
%positive or negative z-direction of its pullback to 2D
%   Input Parameters:
%       - mesh:             The mesh struct describing the surface
%       - v2D:              The (U,V)-coordinates of the vertices in the
%                           plane
%
%   Output Parameters:
%       - pm:               The sign of the orientation of the normal in
%                           the plane

f = mesh.f; v3D = mesh.v; vn = mesh.vn;

assert( ~isempty(vn), ...
    'The mesh provided does not have any prescribed normal vectors!');

%--------------------------------------------------------------------------
%       Calculate the Face Normals From the Vertex Normals
%--------------------------------------------------------------------------

% Simply averages the the three vertex normals of a particular face
% fn = sum(permute(reshape(vn(f',:)', 3, 3, []), [2 1 3]), 1);
% fn = reshape( fn, 3, [] )' / 3;
fn = cat( 3, vn(f(:,1),:), vn(f(:,2),:), vn(f(:,3),:) );
fn = sum(fn, 3) ./ 3 ;

% Make sure that the normal vectors have unit norm
fn = fn ./ repmat( sqrt( sum( fn.^2, 2 ) ), 1, 3 );

%--------------------------------------------------------------------------
%    Calculate the Normal Vectors Defined by the Edges of the Faces
%--------------------------------------------------------------------------

% For the mesh in 3D
e1_3D = v3D(f(:,2), :) - v3D(f(:,1), :);
e2_3D = v3D(f(:,3), :) - v3D(f(:,1), :);

nn_3D = cross( e1_3D, e2_3D, 2 );
nn_3D = nn_3D ./ repmat( sqrt( sum( nn_3D.^2, 2 ) ), 1, 3);

%--------------------------------------------------------------------------
%   Determine the Alignment of the Prescribed Normal and the Calculated
%   Normal. Negative if the orientation is opposite of vertex normals.
%--------------------------------------------------------------------------

pm_3D = sign( dot( fn, nn_3D, 2 ) );

%--------------------------------------------------------------------------
% Now for the mesh in 2D
%--------------------------------------------------------------------------
e1_2D = v2D(f(:,2), :) - v2D(f(:,1), :);
e2_2D = v2D(f(:,3), :) - v2D(f(:,1), :);

e1_2D = [ e1_2D, zeros(size(e1_2D,1), 1) ];
e2_2D = [ e2_2D, zeros(size(e2_2D,1), 1) ];

nn_2D = cross( e1_2D, e2_2D, 2 );
nn_2D = nn_2D ./ repmat( sqrt( sum( nn_2D.^2, 2) ), 1, 3);

% The orientation of the normal vector (i.e., plus/minus z-direction);
% If the faces are counterclockwise, then pm_2D will be all positive
pm_2D = sign( nn_2D(:,3) );

%--------------------------------------------------------------------------
%   Determine the orientation of the prescribed normal in the plane
%--------------------------------------------------------------------------

pm = pm_2D .* pm_3D;


end
