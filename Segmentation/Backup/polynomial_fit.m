function [x_fit,y_fit] = polynomial_fit(x, y)
% Input: x and y vector of equal length 
% Output: two vectors delineating x, and y components of polinomial fit for
% OCT Images

% define variables
xx=[6:504];
yy=[];

% polyfit and filter inputs
fit = polyfit(x,y,3);

% acquire y fitted values
for i=1:length(xx)
    yy(i)=fit(4)+fit(3)*xx(i)+fit(2)*xx(i)^2+fit(1)*xx(i)^3;
end

x_fit = xx;
y_fit = yy;
end

