
function sysBal= balancing(sysPss,T,Ti,dT,gridPoints)
%BALANCING Computes the balanced system in pss form (LPVTools)
%with respect to transformation T.
%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------


sysdata = sysPss.Data;

numParams = size(gridPoints,2);
numGridPointsPerParam = nthroot(size(gridPoints,1),numParams);
numGridPointsPerParamRate = 2;

n = size(sysdata.A(:,:,1),1);
m = size(sysdata.B(:,:,1),2);
l = size(sysdata.C(:,:,1),1);

str = [];
for i=1:2*numParams
   
    if(i > numParams)
        str = [str ',' num2str(numGridPointsPerParamRate)]; 
    else 
        str = [str ',' num2str(numGridPointsPerParam)]; 
    end
end

eval(['ABal = NaN(' num2str(n) ',' num2str(n) str ');']);
eval(['BBal = NaN(' num2str(n) ',' num2str(m) str ');']);
eval(['CBal = NaN(' num2str(l) ',' num2str(n) str ');']);
eval(['DBal = NaN(' num2str(l) ',' num2str(m) str ');']);


paramNames = fieldnames(sysPss.Parameter);
paramNames = [paramNames;cell(numParams,1)];

grids = cell(numParams,1);

for i = 1:numParams
   
    paramNames{numParams+i} = [char(paramNames{i}) 'd'];
    grids{i} = pgrid(paramNames{numParams+i},linspace(sysPss.Parameter.(paramNames{1}).Ratebound(1),sysPss.Parameter.(paramNames{1}).Ratebound(2),numGridPointsPerParamRate),[-inf,inf]);
end

str = [];
for i = 1:numParams
   
    str = [str ',grids{' num2str(i) '}'];
end

eval(['newDomain = rgrid(sysPss.Domain' str ');']);


ratePoints = findAllGridPoints(numParams,numGridPointsPerParamRate,[sysPss.Parameter.(paramNames{1}).Ratebound(1),sysPss.Parameter.(paramNames{1}).Ratebound(2)]);

index = 1;

for inRatePoint = 1:size(ratePoints,1)   
    for inGridPoint=1:size(gridPoints,1)
        
        sum = zeros(n);
        
        for inParam = 1:numParams
            
           sum = sum + Ti{inGridPoint} * dT{inGridPoint}{inParam}*ratePoints(inRatePoint,inParam);
        end
        
        ABal(:,:,index) = Ti{inGridPoint}*sysdata.A(:,:,inGridPoint) * T{inGridPoint} - sum;       % eye(n);%
        BBal(:,:,index) = Ti{inGridPoint}*sysdata.B(:,:,inGridPoint);
        CBal(:,:,index) = sysdata.C(:,:,inGridPoint)*T{inGridPoint};
        DBal(:,:,index) = sysdata.D(:,:,inGridPoint);
        
        index = index+1;
    end
end

systrans=ss(ABal,BBal,CBal,DBal);




% sysBal=pss(systrans,sysPss.Domain);
sysBal=pss(systrans,newDomain);
end