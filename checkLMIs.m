
function [lmisSatisfied, eigenValue] = checkLMIs(sys,gramian,gramianFlag)
%CHECKLMIS Validates an affine Gramian grid based
%
%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------

    n = order(sys);
    paramNames = fieldnames(sys.Parameter);
    numParams = numel(paramNames);
    numGridPoints = numel(sys.Data.A)/n^2;
    
    % Kind of Gramian
    if(strcmp(gramianFlag,'c'))

        calcCtrb = true;
    elseif(strcmp(gramianFlag,'o'))

        calcCtrb = false;
    else
        
        error('Flag unspecified')
    end

    [gramian_, dGramian] = gramianFctHandle(gramian);

    lmisSatisfied = [true; true];
    
    breakFlag = false;
    
    for i=1:numGridPoints    
       for j = 1:2 

           sum = 0;

           for k = 1:numParams

               sum = sum + dGramian{k} * sys.Parameter.(paramNames{k}).RateBound(j);
           end

           if(calcCtrb)

               LMI = -sum...
                   + sys.Data.A(:,:,i) * evalFctHandleAtVec(gramian_,sys.Domain(i))...
                   + evalFctHandleAtVec(gramian_,sys.Domain(i)) * sys.Data.A(:,:,i)'...
                   + sys.Data.B(:,:,i) * sys.Data.B(:,:,i)';
           else

               LMI = sum...
                   + sys.Data.A(:,:,i)' * evalFctHandleAtVec(gramian_,sys.Domain(i))...
                   + evalFctHandleAtVec(gramian_,sys.Domain(i)) * sys.Data.A(:,:,i)...
                   + sys.Data.C(:,:,i)' * sys.Data.C(:,:,i);
           end
           
           if(max(eig(LMI)) > 0)
               
               lmisSatisfied(2) = false;
               eigenValue{2}(i,j) = max(eig(LMI));
               
               breakFlag = true;
               break
           end 
       end
       
       if(breakFlag)
           break
       end
       
       LMI = evalFctHandleAtVec(gramian_,sys.Domain(i));
       
       if(min(eig(LMI)) < 0)
               
           lmisSatisfied(1) = false;
           eigenValue{1}(i,j) = min(eig(LMI));
           break
       end 
       
       
       
    end


    if(all(lmisSatisfied))

       eigenValue = NaN;
    end







end

