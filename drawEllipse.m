
function drawEllipse(W,gramianType)
%DRAWELLIPSE Draws ellipse for a given Gramian and unit energy
%
%---------------------------------------------------------------------------------------------------
% For Paper
% "Grid-Free Constraints for Parameter-Dependent Generalized Gramians via Full Block S-Procedure"
% by Lennart Heeren and Herbert Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Lennart Heeren
%---------------------------------------------------------------------------------------------------

    lineWidth = 1;

    if(size(W,1)>2)
        
        error('W exceed dimension of 2')
    end

    [U,S] = svd(W);

    if(strcmp(gramianType,'o'))
       
        % obsv
        a=sqrt(1/S(1,1)); % horizontal radius
        b=sqrt(1/S(2,2)); % vertical radius
    elseif(strcmp(gramianType,'c'))
        
        % ctrb
        a=sqrt(S(1,1)); % horizontal radius
        b=sqrt(S(2,2)); % vertical radius
    end
    
    phi=-1.05*pi:0.01:1.05*pi;
    coords = [a*cos(phi); b*sin(phi)];    
    
%     coords = linspace(-a,a,1000);
%     coords(2,:) = (b/a) * sqrt(a^2 - coords.^2);    
%     coords = [coords, [coords(1,:); -coords(2,:)]];

    coords = U * coords;
    
    if(strcmp(gramianType,'o'))
        
        for i = 1:size(coords,2); sol(i,1) = coords(:,i)' * U * S * U' * coords(:,i); end
    else
        
        for i = 1:size(coords,2); sol(i,1) = coords(:,i)' * U' / S * U * coords(:,i); end
    end
    
    if(abs(sol-1) > 1e-10)
        
        dlg = warndlg('Ellipse equation is not fullfiled!');
        f = figure(100); plot(sol);
        
        waitfor(dlg)
        close(f)
        wait(.1)
    end
    
    principleAxis = U * diag([a b]);
    
    hold on
    p = plot(coords(1,:),coords(2,:),'LineWidth',lineWidth);
    c = compass(principleAxis(1,:),principleAxis(2,:),'--');
    c(end-1).Color = p(end).Color;
    c(end).Color = p(end).Color;
    set(get(get(c(end-1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(c(end),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold off
    
    axis equal
end