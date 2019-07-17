function [ intersects, allRange, alpha ] = calculate_point_interval( Mj,...
    Bk, eps)
%CALCULATE_POINT_INTERVAL Calculates an optimistic angular interval over
%which a landmark point within the unit disk can be rotated to lie within 
%a target circle
%   Input Parameters:
%       Mj:     The landmark point (a complex value)
%
%       Bk:     The target point (a complex value)
%
%       eps:    The radius of the target circle
%
%   Output Parameters:
%       intersects:     A boolean contingent on whether the target circle
%                       intersects the circle defined by the landmark point
%
%       allRange:       A boolean contingent on whether the relevant
%                       angular interval is [-pi,pi]
%
%       alpha:          The angular interval [alpha1, alpha2]
%                       alpha1 is the angle between the point and the most
%                       CW bounding angle of the target circle.  alpha2 is
%                       the angle between the point and the most CCW
%                       bounding angle of the target circle
%                       -pi <= alpha([1,2]) <= pi

intersects = false;
allRange = false;
alpha = [];

D = abs(Bk); % absolute distance to target point
phi = angle(Bk); % -pi <= phi <= pi

R = abs(Mj); % absolute distance to the landmark point
theta = angle(Mj); % -pi <= theta <= phi

%--------------------------------------------------------------------------
% Determine if the landmark circle is contained within the target circle
%--------------------------------------------------------------------------
if eps > (D + R)
    allRange = true;
    alpha = [-pi, pi];
    return; 
end

%--------------------------------------------------------------------------
% Check if the target circle intersects the landmark circle
% If so proceed with the remainder of the algorithm
%--------------------------------------------------------------------------

b_circ_min = D-eps;
b_circ_max = D+eps;

if (b_circ_min <= R) && (R <= b_circ_max)
    
    intersects = true;
    
    % The angular shift from the location of Bk to the intersection points
    % with the landmark circle
    dPhi = acos( (R^2 + D^2 - eps^2) / (2*R*D) );
    
    alpha = [phi - dPhi - theta, phi + dPhi - theta];
    alpha = wrapToPi(alpha);
   
end


end

