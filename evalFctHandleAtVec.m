
function evalfct = evalFctHandleAtVec(fctHandle,vec)
%EVALFCTHANDLEATVEC Evaluates fctHandle at values listed in vec
%
%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------

    numInputs = length(vec);

    gridStr = num2str(vec(1));
    for k=2:numInputs

        gridStr = [gridStr, ',' num2str(vec(k))];
    end

    eval(['fctEval = fctHandle(' gridStr ');']);
end