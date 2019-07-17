function [ int_list ] = add_interval( list1, list2 )
%ADD_INTERVAL Generates the union of a set of intervals such that
%overlapping intervals are merges together
%   Input Parameters:
%       -list1/list2:   asceding sorted list of 2*n points defining n
%                       intervals

if ~isempty(list1)
    A1 = [list1(1:2:end) list1(2:2:end)];
else
    A1 = [];
end

if ~isempty(list2)
    A2 = [list2(1:2:end) list2(2:2:end)];
else
    A2 = [];
end

A = [A1; A2];

n = size(A,1);
[t,p] = sort(A(:));
z = cumsum(accumarray((1:2*n)', 2*(p<=n)-1));
z1 = [0; z(1:end-1)];

int_list = [t(z1==0 & z>0), t(z1>0 & z==0)]';
int_list = int_list(:);


end

