function [ theta, maxStabs ] = optimal_stabbing_angle( S )
%OPTIMAL_STABBING_ANGLE Uses a stabbing query to determine the optimal
%angle that brings the most landmark points within their various target
%circles
%   Input Parameters:
%       S:      A cell array holding the relevant intervals.  The index of
%               the cell array S{j} corresponds to a particular landmark
%               point Mj.  S{j} is an ordered list of the end points of the
%               intervals for Mj (i.e. size(S{j},1) = number of
%               intervals*2).  The second column of S{j} holds a value +1
%               (-1) denoting each point as the start (end) of an interval
%
%   Output Parameter
%       theta:      The optimal angle returned by the stabbing query
%
%       maxStabs:   The maximum number of intervals stabbed by theta

%--------------------------------------------------------------------------
% Query the total number of interval points
%--------------------------------------------------------------------------
numPoints = 0;
for j=1:length(S)
    numPoints = numPoints + size(S{j},1);
end

%--------------------------------------------------------------------------
% Build the data strucure used for the stabbing query
%--------------------------------------------------------------------------
stabStruct = zeros(numPoints,2);
ind = 1;
for j=1:length(S)
    
    len = size(S{j},1);
    stabStruct(ind:(ind+len-1),:) = S{j};
    ind = ind+len;

end

% Sort the structure first by end points in ascending order and then by
% interval value in descending order (so that start points precede end
% points with the same location)
stabStruct = sortrows(stabStruct, [1, -2]);

% Build the running sum of the interval values
stabStruct = [stabStruct cumsum(stabStruct(:,2))];

%--------------------------------------------------------------------------
% Determine the optimal location
%--------------------------------------------------------------------------

maxStabs = max(stabStruct(:,3));

% If the stabbing problem is poorly defined, there may be multiple
% locations that stab the same number of intervals.  Here we find all such
% 'optimal' locations
maxPnts = find( stabStruct(:,3) == maxStabs );

% Next we optimize the locations to lie in the middle of their appropriate
% interval instead of along an end point and break the degeneracy by
% choosing the point with the largest interval
% if length(maxPnts) > 1
%     warning(['Stabbing problem not well defined. Breaking ',...
%         num2str(length(maxPnts)), '-fold degeneracy']);
% end

maxInt = 0;
theta = 0;

for i = 1:length(maxPnts)
    
    startPnt = stabStruct(maxPnts(i),1);
    endID = find( (stabStruct(:,1) > startPnt) & ...
        (stabStruct(:,2) == -1), 1);
    endPnt = stabStruct(endID);
    
    if (endPnt - startPnt) > maxInt
        maxInt = endPnt - startPnt;
        theta = (endPnt - startPnt)/2 + startPnt;
    end
    
end


end

