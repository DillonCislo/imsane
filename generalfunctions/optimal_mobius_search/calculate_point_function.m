function [ Uz, theta ] = calculate_point_function( z, M, B, eps )
%CALCULATE_POINT_FUNCTION Calculate the optimization function for a single
%Mobius parameter
%   Input Parameters:
%       z:      The Mobius parameter defining a transformation (a complex
%               number)
%
%       M:      The landmark point set
%
%       B:      The target point set
%
%       eps:    The radius of the target circles

% This cell array holds the angular ranges for which transformed landmark
% point intersects a target circle
ang_ranges = cell(size(M));

%--------------------------------------------------------------------------
% Cycle through all landmark points to extract the various angular ranges
%--------------------------------------------------------------------------
for j = 1:length(M)
    
    Mj = M(j);
    
    %----------------------------------------------------------------------
    % Calculate the Mobius transform of the current landmark point
    %----------------------------------------------------------------------
    MMj = (Mj - z) / (1 - conj(z) * Mj);
    
    % Holds the relevant angular intervals for this landmark point
    ang_range = [];
    
    %----------------------------------------------------------------------
    % Cycle through the target points for this particular landmark to
    % extract the relevant angular intervals
    %----------------------------------------------------------------------
    Bj = B{j};
    for k = 1:length(Bj)
        
        % Calculate the angular interval in which the transformed landmark
        % point intersects a target circle
        [intersects, allRange, alpha] = calculate_point_interval(MMj, ...
            Bj(k), eps(j));
        
        % If the angular range is [-pi, pi] set the range and continue to
        % the next landmark point
        if allRange
            ang_range = [-pi, pi]';
            continue;
        % If the angular range is a regular subset of [-pi, pi] add the
        % interval to the set of angular ranges for this landmark point
        elseif intersects
            ang_range = add_interval(ang_range, alpha);
        end
        
    end
    
    % Make sure the orientation of ang_range is vertical
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
[theta, Uz] = optimal_stabbing_angle(ang_ranges);


end

