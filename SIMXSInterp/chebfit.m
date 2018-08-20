function y = chebfit(x_data, f_data, x)
%
% y = chebfit(x_data, f_data, x);
%
% Construct and evaluate a Chebyshev representation of the
% polynomial that interpolates the data points (x_i, f_i):
%
% p = b(1)*T_0(x) + b(1)*T_1(x) + ... + b(n)T_N(x)
%
% where n = N+1, and T_j(x) = cos(j*acos(x)) is the jth Chebyshev
% polynomial.
%
n = length(x_data);
xmax = max(x_data);
xmin = min(x_data);
x_data = (2*log10(x_data) - log10(xmax) - log10(xmin))/(log10(xmax) - log10(xmin));

T = zeros(n, n);
T(:,1) = ones(n,1);
T(:,2) = x_data;
for j = 3:n
T(:,j) = 2*x_data.*T(:,j-1) - T(:,j-2);
end
b = T \ log10(f_data);
x = (2*log10(x) - log10(xmax) - log10(xmin))/(log10(xmax) - log10(xmin));
y = zeros(size(x));
for j = 1:n
y = y + b(j)*cos( (j-1)*acos(x) );
end