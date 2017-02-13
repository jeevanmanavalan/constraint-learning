% Function for calculating the Jacobian of a function with finite central differences.
%
%     J = fd ( f, x )
%
% Computes the Jacobian of the function f with finite central differences.
% 
% Author: Jeevan Manavalan & Matthew Howard
%
% input: 
%     f     - function handle to the function
%     x     - point in input space where Jacobian should be evaluated.
%
% output:
%     J - Jacobian
%
function J = fd ( f, x )

delta=1e-6;
for i=1:length(x)
	dx = zeros(size(x)); dx(i) = delta;
	yp = f(x+dx);
	ym = f(x-dx);
	J(:,i) = ((yp(:) - ym(:))/(2*delta))';
end