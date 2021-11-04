
function converted = uss2LpvToolsType(sys,type,numGridPointsPerParam,maxRate)
%USS2LPVTOOLSTYPE Converts an uss model into pss model lft model (both LPV
%tools types)
%
%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------

    if(exist('numGridPointsPerParam','var'))
        
        thirdArg = true;
    else
        
        thirdArg = false;
    end
    
    if(~exist('maxRate','var'))
        
        maxRate = Inf;
    end


    if(isa(sys,'uss'))

        parameterNames = fieldnames(sys.Uncertainty);
        numParameters = numel(parameterNames);

        rateBounds = cell(numParameters,2);

        for i = 1:numParameters

            rateBounds{i,1} = parameterNames{i};
            rateBounds{i,2} = [-maxRate maxRate];
        end

        lpvToolsLFT = plftss(sys,rateBounds);

        if(strcmp(type,'pss'))

            if(thirdArg)
                converted = lft2grid(lpvToolsLFT,numGridPointsPerParam);
            else 
                converted = lft2grid(lpvToolsLFT,3);
            end
        elseif(strcmp(type,'lft'))

            converted = lpvToolsLFT;
        else

            converted = NaN;    
        end
    else

        converted = NaN;
    end
end