function [ US, theta ] = calculate_bounding_function( S, M, B, eps )
%CALCULATE_BOUNDING_FUNCTION Calculate the bounding function for a square
%subset of the unit disk
%   Input Parameters:
%       S:      The struct containing the geometric information defining
%               the square region
%
%       M:      The landmark point set
%
%       B:      The target point sets


% This cell array holds the angular ranges for which the bounding region of
% the map of each element of M intersects a target circle
ang_ranges = cell(size(M));

%--------------------------------------------------------------------------
% Calculate the circumcircle of the square input region
%--------------------------------------------------------------------------
cR = complex(S.Center(1), S.Center(2));
rR = S.Pos(3) / sqrt(2);

%--------------------------------------------------------------------------
% Cycle through all landmark points to extract the various angular ranges
%--------------------------------------------------------------------------
for j = 1:length(M)

    %----------------------------------------------------------------------
    % Calculate the bounding region
    %----------------------------------------------------------------------
    [rMin, rMax, thetaMin, thetaMax] = calculate_bounding_region( cR, ...
        rR, M(j));
    
    % Holds the end points of the angular intervals for this particular
    % point
    ang_range = [];
    
    %----------------------------------------------------------------------
    % Cycle through the target points for this particular landmark to
    % extract the relevant angular intervals
    %----------------------------------------------------------------------
    Bj = B{j};
    for k = 1:length(Bj)
        
        % Calculate the angular interval in which the bounding wedge
        % intersects the target circle defined by Bj(k)
        [intersects, allRange, alpha] = calculate_wedge_interval(Bj(k), ...
            eps(j), rMin, rMax, thetaMin, thetaMax);
        
        % If the angular interval is [-pi,pi] set the range and continue to
        % the next landmark point
        if allRange
            ang_range = [-pi; pi];
            continue;
        % If the angular range is a regular subset of [-pi, pi] add the
        % interval to the set of angular ranges for this landmark point
        elseif intersects
            ang_range = add_interval(ang_range, alpha);
        end

    end
    
    % Make sure that the orientation of ang_range is vertical
    if size(ang_range,2) > size(ang_range,1)
        ang_range = ang_range';
    end
    
    intVal = zeros(length(ang_range),1);
    intVal(1:2:end) = 1;
    intVal(2:2:end) = -1;
    
    ang_range = [ang_range intVal];
    
    % Add the current angular range to the storage array
    ang_ranges{j} = ang_range;
    
end

%--------------------------------------------------------------------------
% Determine the optimal stabbing angle
%--------------------------------------------------------------------------
[theta, US] = optimal_stabbing_angle(ang_ranges);

end

