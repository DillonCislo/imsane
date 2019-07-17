function [ intersects ] = circle_square_intersect( Cc, R, Cs, W, varargin )
%CIRCLE_SQUARE_INTERSECT Determines whether or not the intersection between
%a circular domain and a square domain is non-empty
%   Input Parameters:
%       Cs:     The center of the square as [x,y]
%       W:      The horizontal width of the square
%       (H:      The vertical height of the square)
%       Cc:     The center of the circle as [x,y]
%       R:      The radius of the circle
%
%   NOTE: IT IS ASSUMED THAT THE SQUARE IS AXIS ALIGNED

H = W;

for i = 1:length(varargin)
    if isa(varargin{i}, 'double')
        continue;
    end
    if ~isempty(regexp(varargin{i}, '^[Rr]ectangle', 'match'))
        H = varargin{i+1};
    end
end

xs = Cs(1); ys = Cs(2);
xc = Cc(1); yc = Cc(2);

intersects = false;

% Calculate the distance between the region centers
D = sqrt( (xs - xc)^2 + (ys - yc)^2 );

% Calculate the angle between the x-axis and the line joining the region
% centers (after shifting the origin to the center of the circle)
u = xs - xc; v = ys - yc;
z = complex(u,v);
theta = angle(z); % Note: the angles here lie between [-pi, pi]

% Find the distance from the center of the square to the edge of the square
% along the line that joins the region centers
if (-3*pi/4 <= theta <= -pi/4) || ...
        (pi/4 <= theta <= 3*pi/4)
    L = abs( (H/2) / sin(theta) );
else
    L = abs( (W/2) / cos(theta) );
end

if (L+R) >= D
    intersects = true;
end


end

