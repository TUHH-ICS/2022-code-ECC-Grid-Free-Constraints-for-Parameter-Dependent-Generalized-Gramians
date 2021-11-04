%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------

clear; clc; close all

%%
fileExtension = '';


%% Parameter
param.m = 1;        % kg
param.d = .75;    % Ns/m
param.k0 = .5;       % N/m
param.range = 30;   % Percent
param.changeRate = 1e-2;   % N/ns

param.c = 1;        % 0 or 1, set equal 1, if all masses are connected to the wall by spring&damper
param.c = max(0,min(1,round(param.c))); % Ensures that c equals 1 or 0


%% Init
numMassesMin = 1;
numMassesMax = 20;
numRhoMin = 1;
numRhoMax = 10;
sysIter = 1;
maxOrder = 1;

maxTime = 3600/4;     % sec

results = cell(numMassesMax,numRhoMax);



for numRho = numRhoMin:numRhoMax % Outer loop: number of parameters
    
    if(numRho > numMassesMax)

        break
    end
    
    for numMasses = numMassesMin:numMassesMax % Inner loop: number of masses
        
        if(numRho>numMasses)
            
            results{numMasses,numRho} = [];
        
            % Display progress
            disp([datestr(now) ': n=' num2str(2*numMasses) ', numRho=' num2str(numRho) ': numRho > numMasses -> skipped' ])
        else
        
            % Setup model
            results{numMasses,numRho}.sysUb = setupStateSpace(numMasses,numRho,param);
            [numMeasOut,numCtrIn] = size(results{numMasses,numRho}.sysUb);
            
            % Compute LFT model
            [results{numMasses,numRho}.sysUbLft, results{numMasses,numRho}.Delta, results{numMasses,numRho}.blkStruct] = lftdata(results{numMasses,numRho}.sysUb);
            results{numMasses,numRho}.sysUbLft = ssbal(results{numMasses,numRho}.sysUbLft);
            
            % Start time measure
            tic
            
            % Compute Gramians
            try
                [results{numMasses,numRho}.gramians.P, results{numMasses,numRho}.gramians.Q, results{numMasses,numRho}.gramians.WcNominal, results{numMasses,numRho}.gramians.iterations,results{numMasses,numRho}.note] = ...
                    findTightestGramianBounds(results{numMasses,numRho}.sysUbLft,results{numMasses,numRho}.blkStruct,results{numMasses,numRho}.sysUb,1,param.changeRate*ones(numRho,1));
            catch ME
               
                if(~isempty(strfind(ME.message,'-LMI infeasible')))
                    
                    disp(ME.message)
                    results{numMasses,numRho} = ME.message;
                else
                    
                    rethrow(ME);
                end
            end
            
            if(isstruct(results{numMasses,numRho})) % If Gramians were computed
               
                numGridPoints = 3; % Grid points per parameter
                results{numMasses,numRho}.sysPss = uss2LpvToolsType(lft(results{numMasses,numRho}.Delta,results{numMasses,numRho}.sysUbLft),'pss',numGridPoints,param.changeRate);
                
                % List all combinations of grid points
                gridPoints = findAllGridPoints(numRho,numGridPoints,[results{numMasses,numRho}.sysUb.Uncertainty.rho1.Range(1),results{numMasses,numRho}.sysUb.Uncertainty.rho1.Range(2)]);
                
                % Compute balancing transformations in grid points
                [T, Ti, dT] = balancingTransformation(results{numMasses,numRho}.gramians.Q,results{numMasses,numRho}.gramians.P,gridPoints);
                results{numMasses,numRho}.sysPssBal = balancing(results{numMasses,numRho}.sysPss,T,Ti,dT,gridPoints);
                results{numMasses,numRho}.sysPssTrunc = modred(results{numMasses,numRho}.sysPssBal,3:2*numMasses,'Truncate');
                
                
%                 [results{numMasses,numRho}.sysB,results{numMasses,numRho}.gramians.hsv] = balanceSystem(results{numMasses,numRho}.sysUb,results{numMasses,numRho}.gramians.P, results{numMasses,numRho}.gramians.Q);
                results{numMasses,numRho}.balancingTime = toc;
                
                success = true;
            else
                
                success = false;
            end
           
            
            % Display progress
            if success
                
                message = ['requires ' num2str(round(results{numMasses,numRho}.balancingTime)) 's'];
            else
                
                message = 'failed';
            end
            
            disp([datestr(now) ': n=' num2str(2*numMasses) ', numRho=' num2str(numRho) ': ' message])
            
            
            %Save results every ten steps
            if(sysIter >=10)
                
                maxOrder = max(maxOrder, 2*numMasses);
                saveOutput(param,results,2*numMassesMin,maxOrder,numRhoMin,numRho,fileExtension);
                
                sysIter=1;
            end
            sysIter = sysIter+1;
            
            
            % If calculation time exceeds the limit, save results, limit the number of masses and proceed
            % with a higher number of parameters
            if(success && results{numMasses,numRho}.balancingTime > maxTime)
                
                maxOrder = max(maxOrder, 2*numMasses);
                
                if(numRho == 1)
                
                    numMassesMax = maxOrder/2;
                    
                    results = results(1:numMassesMax,1:min([numMassesMax,numRhoMax]));   
                end
                
                saveOutput(param,results,2*numMassesMin,maxOrder,numRhoMin,numRho,fileExtension);
                break
            end
            
        end
    end
    
    % Break if time 2x the time limit is exceeded
    if(results{numMasses,numRho}.balancingTime > 2*maxTime)
            
        break
    end
end


%% Save output

function saveOutput(param,results,minOrder,maxOrder,minRho,maxRho,extension)

    try
        fileName = ['tempResults/temp_' date '_example3_nX=' num2str(minOrder) '-' num2str(maxOrder) '_nRho=' num2str(minRho) '-' num2str(maxRho) extension];
        save(fileName,'param','results','-v7.3')        
    catch ME
        
        dlg = warndlg('Failed to save results');
        waitfor(dlg)
        
        
        pause(.5) 
    end
end