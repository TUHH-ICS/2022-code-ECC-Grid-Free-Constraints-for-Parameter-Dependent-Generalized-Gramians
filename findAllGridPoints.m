
function gridPoints = findAllGridPoints(numParameters,numGridPointsPerParameter,parameterRange)
%FINDALLGRIDPOINTS Lists all grid points of a multi dimensional grid.
%
%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------

    if(~exist('parameterRange','var'))
        
        parameterRange = [-1 1];
    end

    gridPointsPerParameter = NaN(numGridPointsPerParameter,numParameters);
    numGridPoints = numGridPointsPerParameter^numParameters;

    gridPoints = NaN(numGridPoints,numParameters);

    x = cell(numParameters,1);
    outputString = '[x{1}';
    inputString = 'gridPointsPerParameter(:,1)';

    % For all grid points
    for i = 1:numParameters

    %     gridPointsPerParameter(:,i) = linspace(sysUb.Uncertainty.(parameterNames{i}).Range(1),sysUb.Uncertainty.(parameterNames{i}).Range(2),numGridPointsPerParameter)';
        gridPointsPerParameter(:,i) = linspace(parameterRange(1),parameterRange(2),numGridPointsPerParameter)';

        outputString = [outputString ',x{' num2str(i+1) '}'];
        inputString = [inputString ',gridPointsPerParameter(:,' num2str(i+1) ')'];
    end

    outputString = [outputString(1:strfind(outputString,['{' num2str(i+1) '}'])-3) ']' ];
    inputString =  inputString(1:strfind( inputString,['(:,' num2str(i+1) ')'])-24);

    eval([outputString '=ndgrid(' inputString ');']);

    for i = 1:numParameters

        numBlocks = numGridPointsPerParameter^max(0,numParameters-2);

        if(numParameters == 1)

            listSize = numGridPointsPerParameter;
        else

            listSize = numGridPointsPerParameter^2;
        end

        for j = 1:numBlocks

            listOfElements = reshape(x{i}(:,:,j),[],1);

            gridPoints((j-1)*listSize+1:j*listSize,i) = listOfElements;        
        end
    end

end