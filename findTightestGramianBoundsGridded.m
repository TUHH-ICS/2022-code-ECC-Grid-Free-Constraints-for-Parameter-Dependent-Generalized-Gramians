
function [P, Q, nominalGramian, it, info] = findTightestGramianBoundsGridded(sys,gramianOrder,slack,nominalGramian)
%FINDTIGHTESTGRAMIANBOUNDSGRIDDED Minimizes Gramians in a D-K iterative scheme and
%with grid based constraints
%
%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------

    numParams = numel(fieldnames(sys.Parameter));
   
    it = 1;
    info = '';
    maxIt = 10;

    J = zeros(10,1);
    reltol = 1e-3;

    % Start with a weighting matrix containing identities
        
    numBases = 1+gramianOrder*numParams;

    if(~exist('nominalGramian','var'))

        basisStruct = ones(numBases);
        basisStruct(1:end-1,1:end-1) = zeros(numBases-1);
        nominalGramian = kron(basisStruct,eye(order(sys)));
    end
    
    if(~exist('slack','var'))

        slack = zeros(4,1);
    end
    
    P = calcGramianBound(sys,'c',nominalGramian,slack([2,4]));
    if(isnan(P))
        
        info = 'Initial ctrb-LMI infeasible';
        
        feasible = calcGramianBound(sys,'c',nominalGramian,slack([2,4]),true);
            
        if(~isnan(mean(feasible,'all')))

           info = [info ' even though it`s proven to be feasible'];
        end

        error(info)
    end

    % Initiallly
    Q = NaN;

    while true

        % Calc ctrb Gramian
        previousQ = Q;
        Q = calcGramianBound(sys,'o',P,slack([1,3]));
        
        if(isnan(Q))

            if(it == 1)
                
                info = 'Initial obsv-LMI infeasible';
            else
                
                info = ['Obsv-LMI infeasible in iteration ' num2str(it)];
            end                
            
            feasible = calcGramianBound(sys,'o',P,slack([1,3]),true);
            
            if(~isnan(mean(feasible,'all')))
                
               info = [info ' even though it`s proven to be feasible'];
            end
           
            if(it == 1)
                
                error(info);
            else
                
                disp(info)
                Q = previousQ;
                it = it-1;
                break
            end
        end
        
        % Compute cost
        J(it) = costFcn(P,Q,sys.Domain,numel(sys.Data.A)/order(sys)^2);

        % Check stopping conditions  
        if(it == 1)
            reldiff = 2*reltol;
        else
            reldiff = abs((J(it)-J(it-1))/J(it) );
        end

        if(it >= maxIt) || (reldiff<reltol)

            if(it >= maxIt)
                
                info = 'Maximum number of iterations exceeded';
            end
            break
        end

        it = it+1;
        
        % Compute obsv Gramian
        previousP = P;
        P = calcGramianBound(sys,'c',Q,slack([2,4]));

        if(isnan(P))

            info = ['Ctrb-LMI infeasible in iteration ' num2str(it)];
            
            feasible = calcGramianBound(sys,'c',Q,slack([2,4]),true);
            
            if(~isnan(mean(feasible,'all')))
                
               info = [info ' even though it`s proven to be feasible'];
            end
            
            disp(info)
            P = previousP;
            it = it-1;
            break
        end
    end
end

function gramian = calcGramianBound(sys,gramianFlag,gramWeight,slack,onlyFeasibility)

% System specs
n = order(sys);
paramNames = fieldnames(sys.Parameter);
numParams = numel(paramNames);
numGridPoints = numel(sys.Data.A)/n^2;

numBases = size(gramWeight,1) / n;
% gramOrder = (numBases - 1) / numParams;


% Kind of Gramian
if(strcmp(gramianFlag,'c'))
    
    calcCtrb = true;
elseif(strcmp(gramianFlag,'o'))
    
    calcCtrb = false;
   
else
    error('Flag unspecified')
end


%% Determin the grammian
% Init yalmip
yalmip('clear');
sdpOptions = sdpsettings('solver', 'mosek','verbose',0);
% sdpOptions = sdpsettings('solver', 'sdpt3','verbose',0);


% Define decision variables

%Gramian Bases
X = sdpvar(n*numBases,n*numBases,'symmetric');
X(1:end-n,1:end-n) = zeros(size(X,1)-n);
    
[X_, dX] = gramianFctHandle(X);


LMIConstraints = [];
% slack = 1e-6;
% slack = 0;

for i=1:numGridPoints    
   for j = 1:2 
       
       
       sum = 0;
       
       for k = 1:numParams
          
           sum = sum + dX{k} * sys.Parameter.(paramNames{k}).RateBound(j);
       end
    
       if(calcCtrb)
          
           constraint = -sum...
                      + sys.Data.A(:,:,i) * evalFctHandleAtVec(X_,sys.Domain(i))...
                      + evalFctHandleAtVec(X_,sys.Domain(i)) * sys.Data.A(:,:,i)'...
                      + sys.Data.B(:,:,i) * sys.Data.B(:,:,i)'...
                      <= -slack(2) * diag(ones(n,1));
       else
           
           constraint = sum...
                      + sys.Data.A(:,:,i)' * evalFctHandleAtVec(X_,sys.Domain(i))...
                      + evalFctHandleAtVec(X_,sys.Domain(i)) * sys.Data.A(:,:,i)...
                      + sys.Data.C(:,:,i)' * sys.Data.C(:,:,i)...
                      <= -slack(2) * diag(ones(n,1));
       end
       
       LMIConstraints = [LMIConstraints, constraint];
   end
   
   LMIConstraints = [LMIConstraints, evalFctHandleAtVec(X_,sys.Domain(i)) >= slack(1) * diag(ones(n,1))];
end



% Solve LMIs

if(exist('onlyFeasibility','var') && onlyFeasibility)
    
    diagnostic = optimize(LMIConstraints,[],sdpOptions);    % Check feasibility
else
    gramWeight_ = gramianFctHandle(gramWeight);
    
    diagnostic = optimize(LMIConstraints,costFcn(X_,gramWeight_,sys.Domain,numGridPoints),sdpOptions);     % Find the tightest bound
end


if(isempty(strfind(diagnostic.info,'Successfully solved (')))
    
    gramian = NaN;
else
    
    gramian = value(X);
end

end


function cost = costFcn(decision,weight,domain,numGridPoints)

    if(isa(decision,'double'))
        
        decision = gramianFctHandle(decision);
    end
    
    if(isa(weight,'double'))
        
        weight = gramianFctHandle(weight);
    end
    

    cost = 0;
    
    for i = 1:numGridPoints
        
        cost = cost + trace(evalFctHandleAtVec(decision,domain(i))*evalFctHandleAtVec(weight,domain(i)));
    end
end
