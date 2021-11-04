
function [u,rho,time] = createExcitation(numIn,numPara,Ts,TFinal,maxValue,maxRate,derivative)
%CREATEEXCITATION Creates normal distributed, uncorrelated signals trhat
%can can be used for the excitation of a LPV model.
%   numIn: Number of system inputs
%   numPara: Number of scheduling parameters
%   Ts: Sampling time
%   TFinal: Length of the sequence
%
%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------

evaluate = false;

    
samplesPerPeriod = TFinal/Ts + 1;
samplesPerPeriodE = round(samplesPerPeriod *1.25);

cutOffPara= min([1, (2*Ts)*maxRate/(2*pi)]);

sequence = idinput([samplesPerPeriodE,numIn+numPara,1],'rgs',[0 1]);
% sequence(:,numIn+1:numIn+numPara) = idinput([samplesPerPeriod,numPara,1],'rgs',[0 cutOffPara]);

sequence(:,1:numIn) = maxValue(1) * sequence(:,1:numIn);
sequence(:,numIn+1:end) = maxValue(2) * sequence(:,numIn+1:end);

for i = numIn+1:numIn+numPara
    
    sequence(:,i) = lowpass(sequence(:,i),cutOffPara,'Steepness',0.99);
end

% Saturate data at three times it's standard deviation
sigma = std(sequence,0,1);
factor = 3;
sequence = min(factor*sigma, max(-factor*sigma, sequence));

% Remove mean
sequence = sequence - mean(sequence,1);

% Scale
for i = 1:size(sequence,2)
    
    if(i>numIn)
        sequence(:,i) = .95 * maxValue(2) * sequence(:,i) ./ max(abs(sequence(:,i)));
    else
        sequence(:,i) = .95 * maxValue(1) * sequence(:,i) ./ max(abs(sequence(:,i)));
    end
end

% Adds derivative
if(exist('derivative','var') &&  derivative)
    
   sequence(2:end-1,numIn+numPara+1:numIn+2*numPara) = (sequence(3:end,numIn+1:numIn+numPara) - sequence(1:end-2,numIn+1:numIn+numPara)) / (2*Ts);    
   sequence = sequence(2:end-1,:);
end

sequence = sequence(end-samplesPerPeriod+1:end,:);

% Separation
u = sequence(:,1:numIn);
rho = sequence(:,numIn+1:end);

% Time vector
time = (0:Ts:(size(sequence,1)-1)*Ts)';



%% Evaluate signal

if(evaluate)
    % Cross correlation
    R = corrcoef(sequence);
    correlationCoeff = max(R-eye(size(R,1)));

    % Mean value
    meanValue = mean(sequence);
    
    % Check for normal distribution
    yesNo = {'Not Normal Distributed', 'Normal Distributed'};
    
    for i = 1:numIn+numPara
        
        normalDistributed(i) = ~kstest(sequence(:,i));
        normalDistributedStr{i} = yesNo(1 + normalDistributed(i));
    end

    plotNoiseSequences(time,sequence,normalDistributedStr,100);
end

end

