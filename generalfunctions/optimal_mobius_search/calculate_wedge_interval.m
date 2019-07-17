function [ intersects, allRange, alpha ] = calculate_wedge_interval( Bk,... 
    eps, rMin, rMax, thetaMin, thetaMax )
%CALCULATE_WEDGE_INTERVAL Calculates an optimistic angular interval over
%which the bounding wedge of a mapped region will intersect with a target
%circle
%   Input Parameters:
%       Bk:         The target point (a complex value)
%
%       eps:        The radius of the target circle
%
%       rMin:       The inner radius of the bounding annulus
%
%       rMax:       The outer radius of the bounding annulus
%
%       thetaMin:   The most CW bounding angle of the wedge
%                   -pi <= thetaMin <= pi
%
%       thetaMax:   The most CCW bounding angle of the wedge
%                   -pi <= thetaMax <= pi
%
%   Output Parameters:
%       intersects:     A boolean contingent on whether the target circle
%                       intersects the bounding annulus
%
%       allRange:       A boolean contingent on whether the relevant
%                       angular interval is [-pi,pi]
%
%       alpha:          The angular interval [alpha1, alpha2]
%                       alpha1 is the angle between the most CCW bounding
%                       angle of the wedge interval and the most CW
%                       bounding angle of the circle.  alpha2 is the angle
%                       between the most CW bounding angle of the wedge and
%                       the most CCW bounding angle of the circle.
%                       -pi <= alpha([1,2]) <= pi


intersects = false;
allRange = false;
alpha = [];

D = abs(Bk); % absolute distance to target point
phi = angle(Bk); % -pi <= phi <= pi

%--------------------------------------------------------------------------
% Determine if one or both of the circles defining the bounding annuli is
% contained in the target circle
%--------------------------------------------------------------------------
in_Target = false;

if eps > (D + rMax)
    in_Target = true;
elseif eps > (D + rMin)
    in_Target = true;
end

if in_Target
    allRange = true;
    alpha = [-pi, pi];
    return;
end

%--------------------------------------------------------------------------
% Check if the target circle intersects either bounding circle
% If so proceed with the remainder of the algorithm
%--------------------------------------------------------------------------

b_circ_max = D+eps;
b_circ_min = D-eps;

if (( (b_circ_min <= rMax) && (rMax <= b_circ_max) ) || ...
        ( (b_circ_min <= rMin) && (rMin <= b_circ_max) ) ) || ...
        ( ( (rMin <= b_circ_min) && (b_circ_min <= rMax) ) && ...
          ( (rMin <= b_circ_max) && (b_circ_max <= rMax) ) ) 
        
    
    intersects = true;
    
    %----------------------------------------------------------------------
    % Check if the angular range of the bounding wedge is [-pi, pi]
    % If so, return the full range
    %----------------------------------------------------------------------
    if (thetaMin == -pi) && (thetaMax == pi)
        allRange = true;
        alpha = [-pi,pi];
        return;
    end
    
    %----------------------------------------------------------------------
    % Extract the bounding interval by finding the tangent lines of the
    % circle originating at the origin
    %----------------------------------------------------------------------
    if D > eps
        
        dPhi = asin(eps/D);
        
        
    %----------------------------------------------------------------------
    % Extract the bounding interval by finding the intersection points with
    % the bounding annulus
    %----------------------------------------------------------------------
    else
        
        % Find intersection points with inner circle (if any)
        dPhiMin = NaN;
        if (b_circ_min <= rMin) && (rMin <= b_circ_max)
            dPhiMin = acos( (rMin^2 + D^2 - eps^2) / (2*rMin*D) );
        end
        
        % Find intersection points with inner circle (if any)
        dPhiMax = NaN;
        if (b_circ_min <= rMax) && (rMax <= b_circ_max)
            dPhiMax = acos( (rMax^2 + D^2 - eps^2) / (2*rMax*D) );
        end
        
        dPhi = max(dPhiMin, dPhiMax);
        
    end
    
    phiMax = phi+dPhi; % The most CCW bounding angle of the circle
    phiMin = phi-dPhi; % The mose CW bounding angle of the circle
    
    alpha = [phiMin - thetaMax, phiMax - thetaMin];
    alpha = wrapToPi(alpha);
    
end

end

