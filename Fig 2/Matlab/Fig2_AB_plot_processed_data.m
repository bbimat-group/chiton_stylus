
%% 
% This script is part of a collection used to visualize data associated
% with a manuscript published in PNAS entitled:
%
% "Persistent polyamorphism in the chiton tooth: 
% from a new biomineral to inks for additive manufacturing."
% 
% by Linus Stegbauer, E. Ercan Alp, Paul J. M. Smeets, Robert Free, Shay G. Wallace, Mark C. Hersam, and Derk Joester* 
%
% DOI:
% 
% This script plots data underlying panels A and B of Figure 2
%
% This script was written by: Maya Kompella and Derk Joester

%% clean up
clear variables
close all
clc

%% Declare constants and flags
expt     = {'SMS Intensity map','SMS Magnetite map'};
id       = {'Fe','Mt'};
fn_inp   = {'Fig2_AB_Femap_processed.csv','Fig2_AB_Mtmap_processed.csv'};
fn_out   = {'Fig2_AB_raw','Fig2_AB_smoothed_contour'};

% set current folder to folder containing the script
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

% give pasth to data and for output relative to current folder
pn_csv = '../csv/';
pn_out = '../';

% flags
save_figure = false;

%% Import data 
for ii = 1:length(fn_inp)
    T{ii} = readtable([pn_csv,fn_inp{ii}]);
end

%% Metadata and formatting:

% unit of imported axes is mm
% iron concentration is normalized to maximum

% define image dimensions
N_rows = 67;
N_cols = 41;

% define coordinate unit
dx = 1000; % [µm]
dy = dx;

% standard deviation of gaussian kernel for smoothing
sigma = 1;

% min and max for color maxi
cmin = [0,0]; 
cmax = [100,100];


%% plot 
close all
clc
for ii = 1:length(expt)
    
    % reshape data to image arrays
    A = reshape(T{ii}{:,1:3},[N_rows,N_cols,3]);
    
    % extract x and y coordinate vectors
    xv = dx*A(1,:,1).'; % [µm]
    xv=xv-min(xv);
    yv = dy*A(:,1,2).'; % [µm]
    yv=yv-min(yv);
    
    % extract intensity data
    I{ii}  = flipud(A(:,:,3));
    
    % normalize to max value and express in %
    In{ii} = I{ii}/max(I{ii},[],'all')*100;
    
    % plot relative intensity as image
%     fig{1,ii} = figure;
% 
%     imagesc(xv([1,end]),yv([1,end]),In{ii});
%     ax = gca;
%     axis xy tight equal
%     ax.TickDir = 'out';
%     caxis(ax,[cmin(ii),cmax(ii)]);
% 
%     cbar = colorbar(ax);
%     cmap = colormap(turbo(cmax(ii)-cmin(ii)+1));
%     cmap(1,:) = [0,0,0];
%     colormap(cmap);
%     cbar.Label.Interpreter = 'latex';
%     cbar.TickDirection = 'out';
%     cbar.Label.FontSize = 14;
%     cbar.Label.String = ['rel. ',id{ii},' conc. [a.u.]'];

 
    % plot relative intensity as filled contour plot
    X = dx*A(:,:,1); %[µm]
    X = X-min(X,[],'all');
    Y = dy*A(:,:,2); %[µm]
    Y = Y-min(Y,[],'all');

    fig{2,ii} = figure;
    
    % smooth with gaussian
    Ins{ii} = imgaussfilt(In{ii},sigma);
    
    % plot filled contour 
    contourf(X,Y,Ins{ii},64,'LineStyle','none');
    axis equal tight
    ax=gca;
    ax.TickDir = 'out';
    cbar2 = colorbar(ax);

    cmap = colormap(turbo(cmax(ii)-cmin(ii)+1));
    cmap(1,:) = [0,0,0];
    colormap(cmap);
    cbar2.Label.Interpreter = 'latex';
    cbar2.TickDirection = 'out';
    cbar2.Label.FontSize = 14;
    cbar2.Label.String = ['rel. ',id{ii},' conc. [a.u.]'];
    caxis(ax,[cmin(ii),cmax(ii)]);

    title({'contour plot','data smoothed w/gaussian kernel',['\sigma = ',num2str(sigma)]}); 
end

% superpose Mt and/or ratio contourlines on normalized Fe intensity map
figure(fig{2,1});
hold on;
contour(X,Y,Ins{2},[15,15],'LineWidth',0.5,'Color','w');
figure(fig{2,2});
hold on;
contour(X,Y,Ins{2},[15,15],'LineWidth',0.5,'Color','w');



if save_figure
    for ii=1:length(expt)
        %saveas(fig{1,ii},[pn_out,fn_out{1},'_',id{ii}],'epsc');
        %disp(['Saved to file:',pn_out,fn_out{1},'_',id{ii}] );
        saveas(fig{2,ii},[pn_out,fn_out{2},'_',id{ii}],'epsc');
        disp(['Saved to file:',pn_out,fn_out{2},'_',id{ii}] );
    end
else
    disp('Figure not saved!');
end