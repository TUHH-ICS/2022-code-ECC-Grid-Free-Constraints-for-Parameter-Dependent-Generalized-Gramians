%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------

clc
clear
close all

%% Setup model

rho = cell(4,1);
for i = 1:numel(rho)

    rho{i} = ureal(['rho' num2str(i)],0,'Range',[-1 1]*.6);
end

GU = ss([-1 .2;0 -1],[1 + rho{1}; 1+rho{2}]*10,[1+rho{3}, 1+rho{4}]*10,0);
param.changeRate = ones(numel(rho),1)*.2;


param.mappedChangeRate = ones(numel(rho),1).* (param.changeRate*2/(GU.Uncertainty.rho1.Range(2)-GU.Uncertainty.rho1.Range(1)));



%% Balance grid based

% Define grid
numGridPoints = 3;
gridPoints = findAllGridPoints(numel(fieldnames(GU.Uncertainty)),numGridPoints,[GU.Uncertainty.rho1.Range(1),GU.Uncertainty.rho1.Range(2)]);
GUgridded = uss2LpvToolsType(GU,'pss',numGridPoints,param.changeRate(1));      

tic % Start time measurement

[WcFbspVaryGridded,WoFbspVaryGridded,~,it] = findTightestGramianBoundsGridded(GUgridded,1,[0;0;1e-5;1e-5]);
GUgriddedCheck = uss2LpvToolsType(GU,'pss',9,param.changeRate(1));


% Validate solution
[lmisWoSatisfied, eigWo] = checkLMIs(GUgriddedCheck,WoFbspVaryGridded,'o');
[lmisWcSatisfied, eigWc] = checkLMIs(GUgriddedCheck,WcFbspVaryGridded,'c');

if(~(all(lmisWcSatisfied) && all(lmisWoSatisfied)))
    
    disp('LMIs not satisfied on a denser grid!!!')
end

balancingTimeVarGridded = toc; %stop time measurment

% Compute balancing transformation
[TVarGridded, TiVarGridded, dTGridded] = balancingTransformation(WoFbspVaryGridded,WcFbspVaryGridded,gridPoints);    

% Balance model
GBVarGridded = balancing(GUgridded,TVarGridded,TiVarGridded,dTGridded,gridPoints);




%% Balance LFT based

% Compute LFT model
[sysGTLft, Delta, blkStruct] = lftdata(GU);

tic % Start time measurement

% Compute parameter dependend gerneralized Gramians
[WcFbspVary,WoFbspVary] = findTightestGramianBounds(sysGTLft,blkStruct,GU,1,param.mappedChangeRate);

balancingTimeVar = toc; %stop time measurment

% Compute balancing transformation
[TVar, TiVar, dT] = balancingTransformation(WoFbspVary,WcFbspVary,gridPoints);

% Balance model
GBVar = balancing(GUgridded,TVar,TiVar,dT,gridPoints);

tic % Start time measurement
% Compute constant gerneralized Gramians
[WcFbspConst,WoFbspConst,~,it] = findTightestGramianBounds(sysGTLft,blkStruct,GU,0);

balancingTimeConst = toc; %stop time measurment

% Balance model
[GB, G] = balanceSystem(GU,WcFbspConst,WoFbspConst);
GBConst = uss2LpvToolsType(GB,'pss',numGridPoints);


%% Truncate models

% Truncate models
GTVarGridded = modred(GBVarGridded,2,'Truncate');
GTVar = modred(GBVar,2,'Truncate');
GTConst = truncateSys(GBConst,1);
    
    
%% Noise response

% Define scheduling trajectory
offset = (GU.Uncertainty.rho1.Range(1)+GU.Uncertainty.rho1.Range(2))/2;
amplitude = (GU.Uncertainty.rho1.Range(2)-GU.Uncertainty.rho1.Range(1))/2;


% Create excitation signal
[uNoise, rhoNoise, tNoise] = createExcitation(1,1,.01,100,[1,(rho{1}.Range(2)-rho{1}.Range(1))/2],param.changeRate(1),true);

ptrajUb.time=tNoise;
ptrajUb.rho1=amplitude*sin(GUgridded.Parameter.rho1.RateBounds(2)/amplitude * tNoise)+offset;
ptrajUb.rho2=amplitude*sin(GUgridded.Parameter.rho1.RateBounds(2)/amplitude * tNoise + .5*pi)+offset;
ptrajUb.rho3=amplitude*sin(GUgridded.Parameter.rho1.RateBounds(2)/amplitude * tNoise + pi)+offset;
ptrajUb.rho4=amplitude*sin(GUgridded.Parameter.rho1.RateBounds(2)/amplitude * tNoise + 1.5*pi)+offset;

ptrajB = ptrajUb;
ptrajB.rho1d= GUgridded.Parameter.rho1.RateBounds(2) * cos(GUgridded.Parameter.rho1.RateBounds(2)/amplitude * tNoise);
ptrajB.rho2d= GUgridded.Parameter.rho1.RateBounds(2) * cos(GUgridded.Parameter.rho1.RateBounds(2)/amplitude * tNoise + .5*pi);
ptrajB.rho3d= GUgridded.Parameter.rho1.RateBounds(2) * cos(GUgridded.Parameter.rho1.RateBounds(2)/amplitude * tNoise + pi);
ptrajB.rho4d= GUgridded.Parameter.rho1.RateBounds(2) * cos(GUgridded.Parameter.rho1.RateBounds(2)/amplitude * tNoise + 1.5*pi);


% Simulate responses
[Y,TUb,X] = lpvlsim(GUgridded,ptrajUb,uNoise,tNoise);

[YBalVar,TBalVar,XBalVar] = lpvlsim(GBVar,ptrajB,uNoise,tNoise);
[YTruncVar,TTruncVar,XTruncVar] = lpvlsim(GTVar,ptrajB,uNoise,tNoise);

[YBalVarGridded,TBalVarGridded,XBalVarGridded] = lpvlsim(GBVarGridded,ptrajB,uNoise,tNoise);
[YTruncVarGridded,TTruncVarGridded,XTruncVarGridded] = lpvlsim(GTVarGridded,ptrajB,uNoise,tNoise);

[YBalConst,TBalConst,XBalConst] = lpvlsim(GBConst,ptrajUb,uNoise,tNoise);
[YTruncConst,TTruncConst] = lpvlsim(GTConst,ptrajUb,uNoise,tNoise);


% Plot response
figure(1)   
subplot(2,1,1)
p = plot(TUb,Y...
        ,TBalVar,YBalVar,TTruncVar,YTruncVar...             
        ,TBalVarGridded,YBalVarGridded,TTruncVarGridded,YTruncVarGridded...
        ,TBalConst,YBalConst,TTruncConst,YTruncConst);
p(1).LineWidth = 2;
p(1).Color = [0.65,0.65,0.65];
p(2).Color = [0.00,0.45,0.74];
p(3).Color = [0.85,0.33,0.10];
p(3).LineStyle = '--';
grid on
title('Noise Response')
legend('unbalanced'...
      ,'balanced (var)','truncated (var)'...           
      ,'balanced (varG)','truncated (varG)'...
      ,'balanced (const)','truncated (const)'...
      ,'Location','eastoutside')

subplot(2,1,2)
plot(TUb,(Y-YBalVar).^2,TUb,(Y-YTruncVar).^2 ...            
    ,TUb,(Y-YBalVarGridded).^2,TUb,(Y-YTruncVarGridded).^2 ...
    ,TUb,(Y-YBalConst).^2,TUb,(Y-YTruncConst).^2)
grid on
title(['Squarred Error:      RMS (const) = ' num2str(round(mean((Y-YTruncConst).^2),2,'significant'),'%10.1e\n')...
       ',      RMS (var) = ' num2str(round(mean((Y-YTruncVar).^2),2,'significant'),'%10.1e\n')...
       ',      RMS (varG) = ' num2str(round(mean((Y-YTruncVarGridded).^2),2,'significant'),'%10.1e\n')])
legend('balanced (var)','truncated (var)'...              
      ,'balanced (varG)','truncated (varG)'...
      ,'balanced (const)','truncated (const)'...
      ,'Location','eastoutside')

