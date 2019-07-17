function [ z, theta ] = optimal_mobius_search( M, B, varargin )
%OPTIMAL_MOBIUS_SEARCH Implementation of the algorithm from 'Conformal
%Surface Alignment Optimal Mobius Search' by Huu Le, Tat-Jun Chin, and
%David Suter (IEEE Conference on Computer Vision and Pattern Recognition
%2016).  Given a set of landmarks and target, the algorithim
%finds an optimal Mobius automorphism of the unit disk that maximizes the
%number of correspondences:
%
%   f(w) = exp(1i*theta) * (w - z) / ( 1 - conj(z) * w )
%
%   Input Arguments:
%       -M:             A vector containing the landmark points
%       -B:             A cell array holding the target points for each
%                       landmark
%
%   Output Arguments
%       - z:              Parameter of optimal Mobius map
%       - theta:          Parameter of optimal Mobius map

%==========================================================================
%                    Pre-Search Data Processing
%==========================================================================

eps = 0.005 .* ones(size(M)); % Radii of the target circles
minSide = 1e-2; % Minimum size of sampling square

dispSearch = true; % Boolean to display search progress

%--------------------------------------------------------------------------
%   Declare Stopping Criteria for Search
%--------------------------------------------------------------------------

MaxIterationNumber = Inf;

% The maximally optimal transformation is one in which each landmark point
% intersects a target circle
Umax = length(M);


for i = 1:length(varargin)
    if isa(varargin{i}, 'double')
        continue;
    end
    if isa(varargin{i}, 'logical')
        continue;
    end
    if ~isempty(regexp(varargin{i}, '^[Rr]adii', 'match'))
        eps = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Mm]in[Ss]ide', 'match'))
        minSide = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Ii]ter[Nn]um', 'match'))
        MaxIterationNumber = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Mm]ax[Mm]atch[Nn]um', 'match'))
        Umax = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Dd]isplay', 'match'))
        dispSearch = varargin{i+1};
    end
end

if isscalar(eps)
    eps = eps .* ones(size(M));
end

assert( isequal(size(eps), size(M)), ...
    ['Radii of target circles must be input either as a scalar ', ...
    'or a vector with the same dimensions as the landmark vector']);

%--------------------------------------------------------------------------
% Construct the tight bounding square of the unit disk
%--------------------------------------------------------------------------
S = struct( 'Pos', [-1 -1 2 2], 'Center', [0,0], ...
    'Priority', Inf);
%--------------------------------------------------------------------------
% Initialize the priority queue and search parameters
%--------------------------------------------------------------------------
Q = S;

% The optimality function and Mobius parameters of the current
% transformation
U = -Inf; z = []; theta = [];

%==========================================================================
%                   RUN OPTIMAL MOBIUS SEARCH!!
%==========================================================================

iter = 0; % clc;
if dispSearch
    fprintf('  Iteration   Optimality  Queue\n');
end

while ( numel(Q) ~= 0 ) && (iter <= MaxIterationNumber)
    
    iter = iter + 1;
    
    %----------------------------------------------------------------------
    % Obtain the square with the highest priority from Q and remove it from
    % the queue
    %----------------------------------------------------------------------
    S = Q(1);
    Q(1) = [];
    
    %----------------------------------------------------------------------
    % Check if the square intersects the unit disk
    %----------------------------------------------------------------------
    intersects = circle_square_intersect( [0,0], 1, S.Center, S.Pos(3) );
    
    if intersects
        
        z0 = complex( S.Center(1), S.Center(2) );
        
        %------------------------------------------------------------------
        % Calculate the optimality function given the parameter z0
        %------------------------------------------------------------------
        [Uz, phi] = calculate_point_function( z0, M, B, eps );
        
        %------------------------------------------------------------------
        % Update values bases on the point optimality function
        %------------------------------------------------------------------
        if Uz == Umax
            U = Uz; z = z0; theta = phi;
            break;
        elseif Uz > U
            U = Uz; z = z0; theta = phi;
        end
        
        %------------------------------------------------------------------
        % Subdivide the current square into four subsquares
        %------------------------------------------------------------------
        SS = divide_square(S);
        
        for i = 1:4
            
            S = SS(i);
            
            if S.Pos(3) >= minSide                
                %----------------------------------------------------------
                % Calculate the optimality function for the square
                %----------------------------------------------------------
                [US, ~] = calculate_bounding_function(S, M, B, eps);
                
                %----------------------------------------------------------
                % Update the queue accordingly
                %----------------------------------------------------------
                if US > U
                    S.Priority = US;
                    Q = add_to_queue(Q,S);
                end
                
            end   
            
        end
        
    end
    
    if dispSearch
        fprintf('%7.0f %15f %8.0f \n', [iter, U, size(Q,2)]);
    end
    
    if iter > MaxIterationNumber
        warning(['Iteration exceeds ', ...
            num2str(MaxIterationNumber), ...
            ' Automatically terminated']);
        break;
    end
    
end



end

