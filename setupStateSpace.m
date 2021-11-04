
function sys = setupStateSpace(numMasses,numRho,param)
%SETUPSTATESPACE Sets up an uncertain state space model of a mass spring
%damper system
%
%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------

    for i = 1:numRho
    
        % Declare scheduling parameters    
        names{i} =  ['rho' num2str(i)];
        rho{i,1} = ureal(names{i}, 0,'Range',[-param.k0*param.range param.k0*param.range]/100);

        % Play with d
%         rho{i,1} = 0;
%         param.d = ureal('rho1',param.d*100,'Percentage',param.range);
    end

    % A
    if(numMasses > 1)
        
        A21(1,1) = -(2*param.k0+rho{1}+rho{mod(2-1,numRho)+1})/param.m;
        A21(1,2) = (param.k0+rho{mod(2-1,numRho)+1})/param.m;

        for i = 2:numMasses-1

            A21(i,i-1) =  (param.k0+rho{mod(i-1,numRho)+1})/param.m;      
            A21(i,i)   = -((2+param.c)*param.k0+param.c*rho{mod(i-1,numRho)+1}+rho{mod(i-1,numRho)+1}+rho{mod(i,numRho)+1})/param.m;
            A21(i,i+1) =  (param.k0+rho{mod(i,numRho)+1})/param.m;

        end

        A21(numMasses,numMasses-1) =  (param.k0+rho{mod(numMasses-1,numRho)+1})/param.m;
        A21(numMasses,numMasses)   = -((1+param.c)*param.k0+param.c*rho{mod(numMasses-1,numRho)+1}+rho{mod(numMasses-1,numRho)+1})/param.m;

        A22 = eye(numMasses)*(-(1+param.c)*param.d/param.m)...
            + blkdiag(0, eye(numMasses-2), 0)*(-param.d/param.m)...
            + [zeros(numMasses-1,1), eye(numMasses-1); zeros(1,numMasses)]*(param.d/param.m)...
            + [zeros(1,numMasses); eye(numMasses-1), zeros(numMasses-1,1)]*(param.d/param.m);
    else
         A21(1,1) = -(param.k0+rho{1}/param.m);
         A22(1,1) = -param.d/param.m;
    end

    A = [zeros(numMasses,numMasses), eye(numMasses); A21, A22];

    % B
    B = zeros(2*numMasses,1);
    B(end) = 1/param.m;

    % C
    C = zeros(1,2*numMasses);
    C(numMasses) = 1;

    % sys
    sys = ss(A,B,C,0);
end

