% Function for calculating the alpha's from the constraint vectors (theta)
%
% A = t2A([x1,x2,x3...xn ])
%
%Author: Jeevan Manavalan
%
% input: 
%     t - constraint vector theta
% output:
%     A - alpha vector
function [A] = t2Alpha(t)

tLength = length(t); dimensionality = tLength+1;
A(1:dimensionality) = [1]; 
if(dimensionality>2)
    for theta = 1:(tLength)
            A(dimensionality) = A(dimensionality)*sin(t(theta));
    end
    cosLimit = 0;
    for pointer = (dimensionality-1):-1:1
        if(pointer)>1
            for theta = 1:(pointer-1)
                A(pointer) = A(pointer)*sin(t(theta));
            end  
        end
        lastTheta = (tLength-cosLimit);
        A(pointer) = A(pointer)*cos(t(lastTheta));
        cosLimit = cosLimit +1;
    end
elseif(dimensionality>1)
    A(dimensionality) = sin(t(1));
    A(dimensionality-1) = cos(t(1));
else
    error('Dimensionality has to be greater than 1')
end


