function [sysB, G, T, Ti] = balanceSystem(sysUb,P,Q)
%BALANCESYSTEM Computes the balancing transformation and the balanced system
%with respect to constant ctrb Gramian P and constant obsv Gramian Q.
%
% (Moore's algorithm, using notation in section 9.5 in
% Robust Systems by Sanchez-Pena and Sznaier)
%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------
    
    [Vc,Sc] = LOCALschursort(P);
    Scrt = sqrt(diag(Sc));
    T1 = lrscale(Vc',1./Scrt,[]); %T1 = inv(Scrt)*Vc'
    T1i = lrscale(Vc,[],Scrt); %T1i = Vc*Scrt
    Wotil = (T1')\(Q/T1);
    [Vo,So] = LOCALschursort(Wotil);
    Sort = diag(So).^(1/4);
    T2 = lrscale(Vo',Sort,[]); %T2 = Scrt*Vo'
    T2i = lrscale(Vo,[],1./Sort); %T2i = Vo*inv(Scrt)

    % Form final transformation
    T = T2*T1;
    Ti = T1i*T2i;

    % Compute balanced system:
    sysB = ss2ss(sysUb,T);

    % Compute Hankel singular values:
    G = sqrt(diag(So));
end


%% schurSort (LPV-Tools)

function [V,S] = LOCALschursort(W)
   
    [V,S] = schur(W);
    [Ssort,idx] = sort(diag(S),'descend');
    S = diag(Ssort);
    V = V(:,idx);  
end