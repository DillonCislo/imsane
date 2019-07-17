function [ Q ] = add_to_queue( Q, S )
%ADD_TO_QUEUE Adds a set of square regions of interest to a priority queue
%such the the output list is sorted by priority
%   Input Parameters:
%       -Q:     The priority queue (a structure array)
%
%       -S:     The list of square (a structure array)

Qfields = fieldnames(Q);
Qcell = struct2cell(Q); szQ = size(Qcell);
if length(szQ) < 3
    szQ(3) = 1;
end
Scell = struct2cell(S); szS = size(Scell);
if length(szS) < 3
    szS(3) = 1;
end

% The outputs of the struct2cell procedure above are 3D arrays
% for an MxN struct array with P fields the size of the converted cell
% array is PxMxN 

% Convert to a matrix and make each field a column
Qcell = reshape(Qcell, szQ(1), [])';
Scell = reshape(Scell, szS(1), [])';

% Concatenate
Qcell = [Qcell; Scell];

% Add column with square size
sideLength = cell(size(Qcell,1),1);
for i = 1:size(Qcell,1)
    
    curPos = Qcell(i,1);
    curPos = curPos{1};
    sideLength{i} = curPos(3);
    
end

Qcell(:,4) = sideLength;

% Sort by Priority
Qcell = sortrows(Qcell, [-3, -4]);

% Remove the side length column
Qcell(:,4) = [];

szQ = [szQ(1) szQ(2) (szQ(3)+szS(3))];

% Put back in original cell array format
Qcell = reshape(Qcell', szQ);

% Convert back to struct array
Q = cell2struct(Qcell, Qfields, 1);



end

