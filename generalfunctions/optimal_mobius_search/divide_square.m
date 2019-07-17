function [ SS ] = divide_square( S )
%DIVIDE_SQUARE Subdivides a given square region into four equal size square
%subregions
%   Detailed explanation goes here

% The side length of the sub squares
L = S.Pos(3) / 2;

Pos = cell(4,1);
Pos{1} = [ S.Pos(1) S.Pos(2) L L ];
Pos{2} = [ S.Pos(1) (S.Pos(2)+L) L L ];
Pos{3} = [ (S.Pos(1)+L) (S.Pos(2)+L) L L ];
Pos{4} = [ (S.Pos(1)+L) S.Pos(2) L L ];

Center = cell(4,1);
Center{1} = S.Center + [-L/2 -L/2];
Center{2} = S.Center + [-L/2 L/2];
Center{3} = S.Center + [L/2 L/2];
Center{4} = S.Center + [L/2 -L/2];

Priority = cell(4,1);
for i = 1:4
    Priority{i} = 0;
end

SS = struct('Pos', Pos, 'Center', Center, 'Priority', Priority);

end

