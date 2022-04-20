
function [T, Ti, dT] = balancingTransformation(QRaw,PRaw,gridPoints)
%BALANCINGTRANSFOIRMATION Computes the balancing transformation 
%based on affine ctrb Gramian P and affine obsv Gramian Q.
%
%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------

numParams = size(gridPoints,2);
numGridPointsPerParam = nthroot(size(gridPoints,1),numParams);

n = size(QRaw,1)/(numParams+1);


strPart = [];
for i=2:numParams
    
   strPart = [strPart ',' num2str(numGridPointsPerParam)]; 
end

if(numParams > 1)
    
    str = [' = cell(' num2str(numGridPointsPerParam) strPart ');'];
else
    str = [' = cell(' num2str(numGridPointsPerParam) ',1);'];
end



%%

[Q,dQ] = gramianFctHandle(QRaw);
[P,dP] = gramianFctHandle(PRaw);

eval(['U' str]);
eval(['S' str]);
eval(['V' str]);

eval(['T' str]);
eval(['Ti' str]);
eval(['dT' str]);


eval(['Rp' str]);
eval(['dRp' str]);

eval(['Rq' str]);
eval(['dRq' str]);


eval(['Zp2' str]);
eval(['dS' str]);
eval(['dV' str]);






for i = 1:numel(Rp)

    
    Qi = evalFctHandleAtVec(Q,gridPoints(i,:));
    Pi = evalFctHandleAtVec(P,gridPoints(i,:));
    
    if i>1
        [T{i}, Ti{i}, Rq{i}, Rp{i}, S{i}, V{i}, U{i}] = computeT(Qi,Pi,S{i-1}, V{i-1});
    else
        [T{i}, Ti{i}, Rq{i}, Rp{i}, S{i}, V{i}, U{i}] = computeT(Qi,Pi);
    end
    
    
    dT{i} = computeTdPert(Rp{i},Rq{i},S{i},V{i},dQ,dP,n,numParams);
%     dT{i} = computeTdCda(Q,P,gridPoints(i,:),delta);

      
end

end

function [T, Ti, Rq, Rp, S, V, U] = computeT(Q,P,previousS,previousV)

    Rq = chol(Q,'upper');
    Rp = chol(P,'lower');

    [U,S,V] = svd(Rq*Rp);
    
    if(exist('previousS','var'))
        [U,S,V] = checkDistance(U,S,V,previousS,previousV);
    end
    
    T = Rp * V / sqrt(S);
    Ti = sqrt(S) \ U' * Rq;
end

function dT = computeTdPert(Rp,Rq,S,V,dQ,dP,n,numParams)

    dTSub = cell(numParams,1);
    dRpSub = cell(numParams,1);
    dRqSub = cell(numParams,1);
    Zp2Sub = cell(numParams,1);
    dSSub = cell(numParams,1);
    dVSub = cell(numParams,1);

    for j = 1:numParams

        dRqSub{j} = sylvester(Rq',Rq,dQ{j});
        dRpSub{j} = sylvester(Rp,Rp',dP{j});

        Zp2Sub{j} = (dRpSub{j}'* Rq'    * Rp * Rq)...
                  + (Rp'    * dRqSub{j}'* Rq * Rp)...
                  + (dRpSub{j}'* Rq'    * Rp * Rq)'...
                  + (Rp'    * dRqSub{j}'* Rq * Rp)';

        dSSub{j} = zeros(n,n);
        dVSub{j} = zeros(n,n);

        for k = 1:n

            dSSub{j}(k,k) = V(:,k)' * Zp2Sub{j} * V(:,k) / (2 * S(k,k));

            sum = 0;
            for m = 1:n

                if(m~=k)

                    sum = sum + (V(:,m)' * Zp2Sub{j} * V(:,k) / (S(k,k)^2 - S(m,m)^2)) * V(:,m);
                end
            end
            dVSub{j}(:,k) = sum;
        end

        dTSub{j} = dRpSub{j} * (V / sqrt(S))...
                 + Rp * ((dVSub{j} / sqrt(S)) - .5 * (V / S^(1.5)) * dSSub{j}); 

    end

    dT = dTSub;
end

function dT = computeTdCda(Q,P,gridPoints,delta)

    numParams = length(gridPoints);
    dT = cell(numParams,1);
    
    for i = 1:numParams
    
        lower = gridPoints;
        lower(i) = lower(i) - delta(i);
        
        upper = gridPoints;
        upper(i) = upper(i) + delta(i);
        
        Qlower = evaluateFHandleAtVec(Q,lower);
        Qupper = evaluateFHandleAtVec(Q,upper);
        
        PLower = evaluateFHandleAtVec(P,lower);
        PUpper = evaluateFHandleAtVec(P,upper);
        
        TLower = computeT(Qlower,PLower);
        TUpper = computeT(Qupper,PUpper);
        
        dT{i} = (TUpper - TLower) / (2*delta(i));
    end
end

function [mappedU,mappedS,mappedV] = checkDistance(currentU,currentS,currentV,previousS,previousV)

    D = NaN(length(currentS));
    A = D;
    
    for i = 1:length(previousS)
        for j = 1:length(currentS)

            z1 = 1/previousS(i,i);
            z2 = 1/currentS(j,j);
            
            D(i,j) = abs((z1-z2)/(z1*z2));
            A(i,j) = (1-abs(previousV(:,i)'*currentV(:,j)));
        end
    end
    
    C = D.*A;
    M = matchpairs(C,1e6);
    
    
    mappedV = currentV(:,M(:,1));
    
    currentS_ = diag(currentS);
    mappedS = diag(currentS_(M(:,1)));
    
    mappedU = currentU(:,M(:,1));
    
    
    % Avoid sign flips in eigenvectors
    for i = 1:length(previousS)
       
        if(previousV(:,i)'*mappedV(:,i) < 0)
            
            mappedV(:,i) = -mappedV(:,i);
            mappedU(:,i) = -mappedU(:,i);
        end
    end
end