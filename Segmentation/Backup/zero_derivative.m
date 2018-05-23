function [x_res, y_res,indexes] = zero_derivative(x, y)
indexes = [];
for i = 1:length(x) 
    for j = i:length(x)
       if (abs(x(i) - x(j)) <= 0)
            indexes = [indexes, x(i)];
       end
    end
end

x_res = x;
y_res = y;
end

