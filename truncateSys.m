
function sysTruncated = truncateSys(sys,reducedSize)
%TRUNCATESYSTEM Truncates pss (LPV Tools) or uss model
%
%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------

    if(isa(sys,'pss'))
        
        originalSize = order(sys);
        elim = (reducedSize+1):originalSize;

        sysTruncated = modred(sys,elim,'Truncate');
        
    elseif(isa(sys,'uss'))
        
        A = sys.A(1:reducedSize,1:reducedSize);
        B = sys.B(1:reducedSize,:);
        C = sys.C(:,1:reducedSize);
        D = sys.D;

        sysTruncated = ss(A,B,C,D);
    else
        
        error('Type not supported')
    end
end