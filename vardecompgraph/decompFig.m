%---------------------------------------------------------------------------%
%---------------------------------------------------------------------------%
% File: figure_10-11-12.m                                                   
%                                                                           
% Description: This code creates figures 10, 11 and 12 in Ajello's          
% ``Financial Intermediation, Investment Dynamics and Business Cycle        
%   Fluctuations.''                                                         
%                                                                           
% Input: A mat file containing the decomposition of endogenous variables
%        in terms of structural shocks.
%        Currently, this code looks two directories up for the file         
%        Ajello_results.mat.                                                
%        This code also depends on Neil Tandon's hatchfill.m to create 
%        speckling, and Philip Graff's CubeHelix for coloring. 
%                                                                           
% Output: Three pdf files, figure10.pdf, figure11.pdf, and figure12.pdf
%         
%                                                                           
% Author: Miguel Acosta                                                     
%---------------------------------------------------------------------------%
%---------------------------------------------------------------------------%
clear all;

%---------------------------------------------------------------------------%
%----------------------- Read and set up data ------------------------------%
%---------------------------------------------------------------------------%
load vardecomp.mat
% Loop for each decomposed variable
for j=1:1
    % List of shock names (make sure these are in the same order as the 'shocks' char array above'
    shock_names = char([       'Stationary Productivity         ';...
                               'Trend Productivity              ';...
                               'Commodity Price                 ';...
                               'Other                           ']);

    

    % List of colors to use, in RGB format 
    cubehelix    = CubeHelix((14*4),2.6,-1.5,2.3,1.5); % Get colors
    colorsetting = fliplr(cubehelix(2:end,:)')';       % Invert and remove black
    colorsetting = colorsetting([2,11,29,2],:); % Select
    %colorsetting = [1 0 0; 0 1 0; 0 0 1; 1 1 0];
    colorsetting(1,:) = [102 102 255]/255;
    
    % Set the dates
    initial_date = 1900;
    firstyear = 1900;
    
    % Get the decomposition data from the Dynare 
    z = [HD' dgdp'];
    
    
    % The rest of this, with slight modifications, is from the Dynare 
    % file 'graph_decomposition.m.' See the Dynare pages for documentation
    % Where I made modifications, 'MA' will be written. 
    comp_nbr = size(z,2)-1;
    gend = size(z,1);
    freq = 1;
    initial_period = 1901;
    x = 1900:2015;

    
    
    z1 = [HD; dgdp];
    xmin = x(1);
    xmax = x(end);
    ix = z1(1:end-1,:) > 0;
    % MA: Changed the axis settings so that there is an extra 5% on either side
    %     of the absolutely tallest bars. 
    ymax = max(sum(z1(1:end-1,:).*ix));
    ix = z1(1:end-1,:) < 0;
    ymin = min(sum(z1(1:end-1,:).*ix));
    pad = (ymax-ymin)*0.05; 
    ymax = ymax + pad;
    ymin = ymin - pad;
    if ymax-ymin < 1e-6
        continue
    end
    f = figure;
    ax=axes('Position',[0.1 0.1 0.6 0.8]);
    plot(ax,x(2:end),z1(end,:),'k-','LineWidth',2)
    axis(ax,[xmin xmax ymin ymax]);
    hold on;
    %text(xmin+(xmax-xmin)/3, ymax+(ymax-ymin)*.02,'Historical Decomposition of Argentine GDP Growth','FontSize', 30) % MA: Plot the title
    box on ; set(gca, 'Layer', 'top'); % MA: Bring axes up
    for i=1:gend
        i_1 = i-1;
        yp = 0;
        ym = 0;
        for k = 1:comp_nbr
            colorsettingCurr = colorsetting(k,:);
            if k == comp_nbr % MA: This is to ensure that the 'initial value' bars are 
                             % always the same colors
                colorsettingCurr = colorsetting(end,:);
            end
            zz = z1(k,i);
            if zz > 0
                a = fill([x(i) x(i) x(i+1) x(i+1)],[yp yp+zz yp+zz yp],colorsettingCurr);
                a = findobj(a, 'type', 'patch'); % MA: This will allow us to add speckling
                yp = yp+zz;
            else
                a = fill([x(i) x(i) x(i+1) x(i+1)],[ym ym+zz ym+zz ym],colorsettingCurr);
                a = findobj(a, 'type', 'patch'); % MA: This will allow us to add speckling
                ym = ym+zz;
            end
            if mod(k,2) == 1 % MA: Every other color is speckled. 
                hatchfill(a,'speckle',45,.25, colorsettingCurr); 
             end
            hold on;
            
        end
    end
    plot(ax,x(2:end),z1(end,:),'k-','LineWidth',4)
        set(gca,'FontSize',17) % Axis font size
    hold off;
    
    axes('Position',[0.75 0.1 0.2 0.8]);
    axis([0 1 0 1]);
    axis off;
    hold on;
    y1 = 0;
    height = 1/comp_nbr;
    labels = char(shock_names,'Initial values');
    
    for i=1:comp_nbr
        colorsettingCurr = colorsetting(i,:);
        if i == comp_nbr
             colorsettingCurr = colorsetting(end,:);
        end
        a = fill([0.1 0.1 0 0],[y1 y1+0.7*height y1+0.7*height y1],colorsettingCurr);
        a = findobj(a, 'type', 'patch'); % MA: This will allow us to add speckling
        if mod(i,2) == 1 % MA: Every other color is speckled. 
            hatchfill(a,'speckle',45,0.25, colorsettingCurr);
        end

        text(0.2,y1+0.3*height,labels(i,:),'Interpreter','tex', 'FontSize', 15);
        hold on
        y1 = y1 + height;
    end
    hold off;
    % MA: Print the pdf
    set(gcf, 'PaperSize', [15 8.5], 'PaperPositionMode', 'manual', ...
             'PaperPosition', [0 0  15 8.5]);
    print(f,'histDecomp','-dpdf');
end
