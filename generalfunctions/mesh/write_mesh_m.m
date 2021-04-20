function write_mesh_m( fullfilename, tr, varargin )
%WRITE_MESH_M(fullfilename, tr, varargin)
%   Write mesh to PLY file
%
% Parameters
% ----------
% fullfilename
% tr
% varargin
% 
% See also
% --------
% plywrite
%

V3D = tr.Points;
F = tr.ConnectivityList;

V2D = zeros(size(V3D, 1),2);
rgb = zeros(size(V3D, 1),3);
VID = 1:size(V3D,1);
C = [];

for i = 1:length(varargin)
    if isa(varargin{i},'double')
        continue;
    end
    if ~isempty(regexp(varargin{i}, '^[Vv]2[Dd]', 'match'))
        V2D = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Rr][Gg][Bb]', 'match'))
        rgb = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Vv][Ii][Dd]', 'match'))
        VID = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Cc]orners', 'match'))
        C = varargin{i+1};
    end
end

if(exist('fullfilename', 'var') == 0)
    [filename, filefolder] = uiputfile('*.m', 'Write m-file');
    fullfilename = [filefolder filename];
end
[filefolder, filename] = fileparts( fullfilename );

fid = fopen(fullfilename, 'w');

edge_list = edges(tr);

ES = freeBoundary(tr);

if ~isempty(C)
    for i = 1:length(C)
        
        new_edges = edge_list( any(edge_list == C(i), 2), : );
        
        ES = [ES; ...
            new_edges(~ismember(new_edges, ES, 'rows'), :) ];
        
    end
end

write_vertices(fid, V3D, V2D, VID, rgb);

write_faces(fid, F);

write_edges(fid, ES);


fclose(fid);


end

function write_vertices(fid, V3D, V2D, VID, rgb)

for i = 1:size(V3D,1)
    
    formatSpec = 'Vertex %d  %0.4f %0.4f %0.4f {uv=(%0.6f %0.6f) rgb=(%0.6f %0.6f %0.6f)}\n';
    fprintf(fid, formatSpec,...
        VID(i), V3D(i,1), V3D(i,2), V3D(i,3), ...
        V2D(i,1), V2D(i,2), ...
        rgb(i,1), rgb(i,2), rgb(i,3) );
    
end

end

function write_faces(fid, F)

for i = 1:size(F,1)
    
   formatSpec = 'Face %d  %d %d %d\n';
   fprintf(fid, formatSpec, ...
       i, F(i,1), F(i,2), F(i,3));
end
    
end

function write_edges(fid, ES)

for i = 1:size(ES,1)
    
    formatSpec = 'Edge %d %d {sharp}\n';
    fprintf(fid, formatSpec, ...
        ES(i,1), ES(i,2) );
    
end

end



