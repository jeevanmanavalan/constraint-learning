% Function for learning the nullspace projection matrix N(q) from observed actions
% 
% Author: Jeevan Manavalan & Matthew Howard    
% 
%If using this code of any purpose, cite the following paper: TBA
%
% input
%     controlPTrain  : observed action training data of the form controlPTest = N(q) * F(q)
%     controlPTest   : observed action test data of the form controlPTest = N(q) * F(q)
%     constraint     : true constraint for evaluation purposes
%     dimensionality : dimensionality of the data controlPTest
% output
%     optimal: N(q) Matrix
function [optimal, results] = gradientDescent(controlPTrain, controlPTest, constraint,dimensionality)

A =@(t)t2Alpha(t);
N =@(t)(eye(dimensionality)-pinv(A(t))*A(t));

startingEstimates = 20;%no of starting points for initial estimated t's
t = constraint;%real constraint(used only for validation purposes)
P = controlPTrain;
PTest = controlPTest;
U=N(t)*P;
UTest = N(t)*PTest;

% function handles for computing error
E =@(t)trace(U'*A(t)'*A(t)*U);
ETest = @(t)trace(UTest'*A(t)'*A(t)*UTest);
% compute error with true data
disp(['True t: ' num2str(E(t)) ' ']);
disp(t);

arrayOfInitialEstimates = [];
arrayListTempEstimates = [];
arrayOfEstimates = [];
for j=1:startingEstimates
    tp = rand(dimensionality-1, 1) *pi;
    arrayOfInitialEstimates = [arrayOfInitialEstimates tp];
    arrayListTempEstimates = [arrayListTempEstimates tp];
    
    count =0;
    prevTP = tp;
    stuckCounter = 0; %keeps track of how many times the error does not increase or decrease in succession
    countLimit = 100;%limit for the number of iterations (how many steps while moving towards the negative of the gradient)
    while E(tp)<=E(prevTP) && (count < countLimit)
        prevTP = tp;
        tp = tp - 1e-3*fd(E,tp)'; %Moving a step towards negative of the gradient(increase precision by decreasing step size 1e-3)
        arrayListTempEstimates = [arrayListTempEstimates tp];
        count = count +1;
        E(tp);
        if(E(tp)==E(prevTP))
            stuckCounter = stuckCounter+1;
            if(stuckCounter==10)
                disp('stuck');
                break;
            end
        else
            stuckCounter = 0;
        end
        if(mod(count,5000)==0)
            disp(['current iteration count for estimate of t ' num2str(count) ' ']);
        end
    end
    arrayOfOngoingEstimates(j).list= arrayListTempEstimates;
    arrayListTempEstimates = [];
    arrayOfEstimates = [arrayOfEstimates tp];
    Up=N(tp)*P;
end
        
bestTP = arrayOfEstimates(:,1);
for h=2:startingEstimates
    currentTP = arrayOfEstimates(:,h);
    if(E(currentTP)<E(bestTP))
        bestTP = currentTP;
        %disp(['new best TP at pos ' num2str(h)]);
        %disp(bestTP);
    end
end
bestTP


results.errorTrain = E(bestTP);%Error on Test Data
results.ppeTrain = get_ppe (U, N(bestTP), P) ;%PPE on Training Data
results.poeTrain = get_poe (U, N(bestTP), P) ;%POE on Training Data 
results.errorTest = ETest(bestTP);%Error on Training Data
results.ppeTest = get_ppe (UTest, N(bestTP), PTest) ;%PPE on Test Data
results.poeTest = get_poe (UTest, N(bestTP), PTest) ;%POE on Test Data   

results.estimatedList.Initial = arrayOfInitialEstimates;%Randomly initialised theta's 
results.estimatedList.Full = arrayOfOngoingEstimates;%Randomly initialised theta's moving towards negative of the gradient
results.estimatedList.Final = arrayOfEstimates;%Final theta's at the end

results.EstimatedT = bestTP;%Estimate of Theta
results.EstimatedA = t2Alpha(bestTP);%Converting Theta to Alpha

optimal.EstimatedP = N(bestTP);%Converting Theta to Nullspace projection matrix N
end
