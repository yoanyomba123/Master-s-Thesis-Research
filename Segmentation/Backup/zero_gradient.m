function x_fin = zero_gradient(x)
% Function : Removes the portions of a line in the x direction that occur multiple times meaning a vertical increas and
% horizontal gradient of 0
% Inputs: A vector;
% Outputs, A vector with limited repeat


% Unique elements
[uv, ~, id] = unique(x);

% How many of each
n = histcounts(id);

% Keep ones with occurences less than 10
x_temp = sort(x, 'ascend');
sec_min = x_temp(end-1);
x_fin = x(ismember(x, uv(n < sec_min)));
end

