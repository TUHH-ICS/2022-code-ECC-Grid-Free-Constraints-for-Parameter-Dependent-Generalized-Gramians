
function [Q, dQ] = gramianFctHandle(QRaw)
%GRAMIANFCTHANDLE Creates fctHandle allows the evaluation of affine Q and its
%derivative at suitable grid points
%
%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------

    if(all(isfinite(value(QRaw)),'all'))
        
        for i = 1:size(QRaw,1)-1
            
           if(~all(QRaw(1:i+1,1:i+1) == 0,'all'))
               
               break
           end
        end
        
        n = max([2,size(QRaw,1)-i]);
    else
    
        n = size(QRaw,1)-sqrt(length(QRaw(isfinite(value(QRaw)))));
    end
    
    
    
    numParams = size(QRaw,1)/n -1;
    
    if(numParams == 0)
        
        Q = QRaw;
        dQ = zeros(size(QRaw,1),size(QRaw,2));
    else
        
        dQ = cell(numParams,1);    
        names = cell(numParams,1);    
        str1=[]; str2=[];

        for i = 1:numParams

            dQ{i} = QRaw(n*(i-1)+1:n*i,end-n+1:end) + QRaw(n*(i-1)+1:n*i,end-n+1:end)';

            if(i>1)

                names{i} = ['rho' num2str(i)];
                str1 = [str1 ',' names{i}];    
                str2 = [str2 '+ ' names{i} ' * (QRaw(' num2str(n*(i-1)+1) ':' num2str(n*i) ',end-n+1:end) + transpose(QRaw(' num2str(n*(i-1)+1) ':' num2str(n*i) ',end-n+1:end)))'];
            end
        end

        eval(['Q = @(rho1' str1 ') QRaw(end-n+1:end,end-n+1:end) + rho1 * (QRaw(1:n,end-n+1:end) + transpose(QRaw(1:n,end-n+1:end)))' str2 ';']);
    end
end