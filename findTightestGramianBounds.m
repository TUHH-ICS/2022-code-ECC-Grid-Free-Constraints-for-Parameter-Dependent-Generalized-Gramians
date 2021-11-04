
function [P, Q, nominalGramian, it, info] = findTightestGramianBounds(sysUbLft,blkStructPlant,sysUb,gramianOrder,maxRate,nominalGramian)
%FINDTIGHTESTGRAMIANBOUNDS Minimizes Gramians in a D-K iterative scheme and
%with gridless constraints
%
%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------

    % Init
    numParams = numel(blkStructPlant);
    
    if(~exist('maxRate','var') || length(maxRate) ~= numParams)
        
        maxRate = Inf;
    end     

    it = 1;
    info = '';
    maxIt = 10;

    J = zeros(maxIt,1);
    reltol = 1e-3;
    
    Q = NaN;

    
    % Start with a weighting matrix containing identities
    if(gramianOrder == 0 || ~isfinite(max(maxRate)))

        if(~exist('nominalGramian','var'))
            nominalGramian = eye(order(sysUb));
        end
        numBases = 1;
    else

        numBases = 1+gramianOrder*numParams;

        if(~exist('nominalGramian','var'))
            basisStruct = ones(numBases);
            basisStruct(1:end-1,1:end-1) = zeros(numBases-1);

            nominalGramian = kron(basisStruct,eye(order(sysUb)));
        end
    end    

    % Compute initial controllability Gramian
    P = calcGramianBoundFBSP(sysUbLft, blkStructPlant, sysUb.Uncertainty,'c',nominalGramian,maxRate);
    if(isnan(P))
        
        info = 'Initial ctrb-LMI infeasible';
        
        feasible = calcGramianBoundFBSP(sysUbLft, blkStructPlant, sysUb.Uncertainty,'c',nominalGramian,maxRate,true);
            
        if(~isnan(mean(feasible,'all')))

           info = [info ' even though it`s proven to be feasible (numerical issue)'];
        end

        error(info)
    end

    % Begin D-K iteration
    while true

        % Calc ctrb Gramian
        previousQ = Q;
        Q = calcGramianBoundFBSP(sysUbLft, blkStructPlant, sysUb.Uncertainty,'o',P,maxRate);
        
        if(isnan(Q))

            if(it == 1)
                
                info = 'Initial obsv-LMI infeasible';
            else
                
                info = ['Obsv-LMI infeasible in iteration ' num2str(it)];
            end                
            
            feasible = calcGramianBoundFBSP(sysUbLft, blkStructPlant, sysUb.Uncertainty,'o',P,maxRate,true);
            
            if(~isnan(mean(feasible,'all')))
                
               info = [info ' even though it`s proven to be feasible (numerical issue)'];
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
        J(it) = costFcn(P,Q,sysUb.Uncertainty);

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
        P = calcGramianBoundFBSP(sysUbLft, blkStructPlant, sysUb.Uncertainty,'c',Q,maxRate);

        if(isnan(P))

            info = ['Ctrb-LMI infeasible in iteration ' num2str(it)];
            
            feasible = calcGramianBoundFBSP(sysUbLft, blkStructPlant, sysUb.Uncertainty,'c',Q,maxRate,true);
            
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
