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

rho1 = ureal('rho1',0,'Range',[-1 1]*.6);
GU = ss([-1 .2;0 -1],[1 + rho1; 1-rho1]*10,[1+rho1, 1-rho1]*10,0);
param.changeRate = .2;
param.mappedChangeRate = param.changeRate*2/(GU.Uncertainty.rho1.Range(2)-GU.Uncertainty.rho1.Range(1));



%% Compute grid based Gramians
gridPoints = linspace(GU.Uncertainty.rho1.Range(1),GU.Uncertainty.rho1.Range(2),3)';
numGridPoints = size(gridPoints,1);

GUgridded = uss2LpvToolsType(GU,'pss',numGridPoints,param.changeRate);      

[WcFbspVaryGridded,WoFbspVaryGridded,~,it] = findTightestGramianBoundsGridded(GUgridded,1,[0;0;1e-6;1e-6]);
GUgriddedCheck = uss2LpvToolsType(GU,'pss',3,param.changeRate);


% Validate solution
[lmisWoSatisfied, eigWo] = checkLMIs(GUgriddedCheck,WoFbspVaryGridded,'o');
[lmisWcSatisfied, eigWc] = checkLMIs(GUgriddedCheck,WcFbspVaryGridded,'c');

if(~(all(lmisWcSatisfied) && all(lmisWoSatisfied)))
    
    disp('LMIs not satisfied on a denser grid!!!')
end



%% Plot unbalanced Gramians

% Compute LFT model
[GULft, Delta, blkStruct] = lftdata(GU);


% LTI 
plotLtiGramiansInGridPoints(GU,gridPoints);


% Compute parameter dependend gerneralized Gramians
[WcFbspVary,WoFbspVary] = findTightestGramianBounds(GULft,blkStruct,GU,1,param.mappedChangeRate);

WoFbspVaryRho = gramianFctHandle(WoFbspVary);
WcFbspVaryRho = gramianFctHandle(WcFbspVary);

WoFbspVaryGriddedRho = gramianFctHandle(WoFbspVaryGridded);
WcFbspVaryGriddedRho = gramianFctHandle(WcFbspVaryGridded);


% Compute constant gerneralized Gramians
[WcFbspConst,WoFbspConst,~,it] = findTightestGramianBounds(GULft,blkStruct,GU,0);


% Plot generalized Gramians
for i = 1:numGridPoints

    subplot(2,numGridPoints,i)

    drawEllipse(WoFbspConst,'o');
    drawEllipse(WoFbspVaryRho(gridPoints(i)),'o')

    if(i==1)
        ylabel({'observability','Gramians'},'fontweight','bold')
    end
    title(['grid point ' num2str(i)])
    axis([-1 1 -1 1] * .3)

    subplot(2,numGridPoints,i+numGridPoints)

    drawEllipse(WcFbspConst,'c');
    drawEllipse(WcFbspVaryRho(gridPoints(i)),'c')
    if(i==1)
        ylabel({'controllability','Gramians'},'fontweight','bold')
    end

    axis([-1 1 -1 1] * 14)
end

% sgtitle('unbalanced')



%% Plot balanced Gramians

%Compute balancing transformations
% [TVarGridded, TiVarGridded] = balancingTransformation(WoFbspVaryGridded,WcFbspVaryGridded,linspace(GT0.Uncertainty.rho1.Range(1),GT0.Uncertainty.rho1.Range(2),numGridPoints)');    
[TVar, TiVar] = balancingTransformation(WoFbspVary,WcFbspVary,linspace(GU.Uncertainty.rho1.Range(1),GU.Uncertainty.rho1.Range(2),numGridPoints)');    
[GB, G, TUb, Ti] = balanceSystem(GU,WcFbspConst,WoFbspConst);

figure(2) % Balanced

% Balance GRamians

% obsvDiagVarGridded = cell(numGridPoints,1);
% ctrbDiagVarGridded = cell(numGridPoints,1);

obsvDiagVar = cell(numGridPoints,1);
ctrbDiagVar = cell(numGridPoints,1);

obsvDiagConst = Ti'*WoFbspConst*Ti;
ctrbDiagConst = TUb*WcFbspConst*TUb';

for i = 1:numGridPoints
    subplot(2,numGridPoints,i)
    warning off
    drawEllipse(zeros(2),'o')
    warning on
    hold on


%     obsvDiagVarGridded{i} = TVarGridded{i}'*WoFbspVaryGriddedRho(gridPoints(i))*TVarGridded{i};
    obsvDiagVar{i} = TVar{i}'*WoFbspVaryRho(gridPoints(i))*TVar{i};

    % LPV varying balanced

%     drawEllipse(obsvDiagVarGridded{i},'o')

    % LPV constant balanced
    drawEllipse(obsvDiagConst,'o')
    drawEllipse(obsvDiagVar{i},'o')
    hold off
    title(['grid point ' num2str(i)])
    title(['ratioConst=' num2str(round(obsvDiagConst(2,2)/obsvDiagConst(1,1)*100)/100)...
           ',  ratioVar=' num2str(round(obsvDiagVar{i}(2,2)/obsvDiagVar{i}(1,1)*100)/100)])
%                ',  ratioVar=' num2str(round(obsvDiagVarGridded{i}(2,2)/obsvDiagVarGridded{i}(1,1)*100)/100)])

    if(i==1)
        ylabel({'observability','Gramians'},'fontweight','bold')
    end  
    axis([-1 1 -1 1] * .35)

    subplot(2,numGridPoints,i+numGridPoints)
    warning off
    drawEllipse(zeros(2),'c')
    warning on
    hold on 

    ctrbDiagVar{i} = TiVar{i}*WcFbspVaryRho(gridPoints(i))*TiVar{i}';
%     ctrbDiagVarGridded{i} = TiVarGridded{i}*WcFbspVaryGriddedRho(gridPoints(i))*TiVarGridded{i}';


%     drawEllipse(ctrbDiagVarGridded{i},'c')
    drawEllipse(ctrbDiagConst,'c')
    drawEllipse(ctrbDiagVar{i},'c')
    hold off

    title(['ratioConst=' num2str(round(ctrbDiagConst(2,2)/ctrbDiagConst(1,1)*100)/100)...
           ',  ratioVar=' num2str(round(ctrbDiagVar{i}(2,2)/ctrbDiagVar{i}(1,1)*100)/100)])
%                ',  ratioVar=' num2str(round(ctrbDiagVarGridded{i}(2,2)/ctrbDiagVarGridded{i}(1,1)*100)/100)])
    if(i==1)
        ylabel({'controllability','Gramians'},'fontweight','bold')
    end
    axis([-1 1 -1 1] * 15)
end

%     sgtitle('balanced')

%% Balance and truncate models

% Define grid
numGridPoints = 10;%3;
gridPoints = findAllGridPoints(numel(GU.Uncertainty),numGridPoints,[GU.Uncertainty.rho1.Range(1),GU.Uncertainty.rho1.Range(2)]);
GUgridded = uss2LpvToolsType(GU,'pss',numGridPoints,param.changeRate);      


% Compute balancing transformations
[TVar, TiVar, dT] = balancingTransformation(WoFbspVary,WcFbspVary,gridPoints);
[TVarGridded, TiVarGridded, dTGridded] = balancingTransformation(WoFbspVaryGridded,WcFbspVaryGridded,gridPoints);    


% Balance models
GBVarGridded = balancing(GUgridded,TVarGridded,TiVarGridded,dTGridded,gridPoints);
GBVar = balancing(GUgridded,TVar,TiVar,dT,gridPoints);
GBConst = uss2LpvToolsType(GB,'pss',numGridPoints);


% Truncate models
GTVarGridded = modred(GBVarGridded,2,'Truncate');
GTVar = modred(GBVar,2,'Truncate');
GTConst = truncateSys(GBConst,1);
    
    
    
%% Noise response

% Define scheduling trajectory
offset = (GU.Uncertainty.rho1.Range(1)+GU.Uncertainty.rho1.Range(2))/2;
amplitude = (GU.Uncertainty.rho1.Range(2)-GU.Uncertainty.rho1.Range(1))/2;


% Create excitation signal
[uNoise, rhoNoise, tNoise] = createExcitation(1,1,.01,100,[1,(rho1.Range(2)-rho1.Range(1))/2],param.changeRate,true);

ptrajUb.time=tNoise;
ptrajUb.rho1=amplitude*sin(GUgridded.Parameter.rho1.RateBounds(2)/amplitude * tNoise)+offset;

ptrajB = ptrajUb;
ptrajB.rho1d= GUgridded.Parameter.rho1.RateBounds(2) * cos(GUgridded.Parameter.rho1.RateBounds(2)/amplitude * tNoise);


% Simulate responses
[Y,TUb,X] = lpvlsim(GUgridded,ptrajUb,uNoise,tNoise);

[YBalVar,TBalVar,XBalVar] = lpvlsim(GBVar,ptrajB,uNoise,tNoise);
[YTruncVar,TTruncVar,XTruncVar] = lpvlsim(GTVar,ptrajB,uNoise,tNoise);

[YBalVarGridded,TBalVarGridded,XBalVarGridded] = lpvlsim(GBVarGridded,ptrajB,uNoise,tNoise);
[YTruncVarGridded,TTruncVarGridded,XTruncVarGridded] = lpvlsim(GTVarGridded,ptrajB,uNoise,tNoise);

[YBalConst,TBalConst,XBalConst] = lpvlsim(GBConst,ptrajUb,uNoise,tNoise);
[YTruncConst,TTruncConst] = lpvlsim(GTConst,ptrajUb,uNoise,tNoise);


% Plot response
figure(5)   
subplot(3,1,1)
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

subplot(3,1,2)
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

subplot(3,1,3)
plot(ptrajB.time,ptrajB.rho1,ptrajB.time,ptrajB.rho1d)
grid on
title('Parameter Trajectory')
legend('rho1                   .','rho1d','Location','eastoutside')
                
                
%% ##################################################################### %%                


function plotLtiGramiansInGridPoints(sys,gridPoints)

    numGridPoints = length(gridPoints);

    for i = 1:numGridPoints
    
        
        rates = linspace(-.01,.01,2);
        
        
        if(isa(sys,'pss'))
            if(i==1)
                    figure()    % LTI
            end
            
            for j=1:length(rates)
                
                G{i,j} = lpvsubs(sys,{'rho1';'rho1d'},[gridPoints(i);rates(j)]);

                Wo{i,j} = gram(G{i,j},'o');
                Wc{i,j} = gram(G{i,j},'c');

                subplot(2,numGridPoints,i)
                drawEllipse(Wo{i,j},'o');

                subplot(2,numGridPoints,i+numGridPoints)
                drawEllipse(Wc{i,j},'c');
            end
            
        else
        
            G{i} = usubs(sys,'rho1',gridPoints(i));

            Wo{i} = gram(G{i},'o');
            Wc{i} = gram(G{i},'c');

            if(i==1)
                figure()    % LTI
            end

            subplot(2,numGridPoints,i)
            drawEllipse(Wo{i},'o');

            subplot(2,numGridPoints,i+numGridPoints)
            drawEllipse(Wc{i},'c');
        end
    end

end
