
function cost = costFcn(decision,weight,uncertainty)
%COSTFCN Computes cost of a Gramian solutioon with respect to a given
%weight.
%
%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------

    if(isa(decision,'double') || isa(decision,'sdpvar'))
        
        decision = gramianFctHandle(decision);
    end
    
    if(isa(decision,'function_handle'))
    
        if(isa(weight,'double'))

            weight = gramianFctHandle(weight);
        end

        %% ########################### %%
        numGridPointsPerParam = 3;
        %% ########################### %%

        gridPoints = findAllGridPoints(numel(fieldnames(uncertainty))...
                                      ,numGridPointsPerParam...
                                      ,[uncertainty.rho1.Range(1),uncertainty.rho1.Range(2)]);
        numGridPoints = size(gridPoints,1);

        cost = 0;

        for i = 1:numGridPoints

            cost = cost + trace(evalFctHandleAtVec(decision,gridPoints(i,:))*evalFctHandleAtVec(weight,gridPoints(i,:)));
        end
    else
        
        cost = trace(decision*value(weight));
    end
end