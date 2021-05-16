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
% This script plots data underlying panels B-H of Figure 3
%
% This script was written by: Maya Kompella and Derk Joester

%% clean up
close all
clear variables
clc


%% Declare constants

fn_inp = {'Fig3_B_processed.csv','Fig3_C_processed.csv','Fig3_D_processed.csv','Fig3_F_processed.csv','Fig3_G_processed.csv','Fig3_H_processed.csv'};
fn_out = 'Fig3B-H';

% flags for output
save_figure = true; 

% set current folder to folder containing the script
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

% give pasth to data and for output relative to current folder
pn_csv = '../csv/';
pn_out = '../';

%% Import data from text files
for ii=1:length(fn_inp)
    opts = detectImportOptions([pn_csv,fn_inp{ii}]);
    opts.VariableTypes = repmat("double",1,length(opts.VariableNames));
    opts = setvaropts(opts, opts.VariableNames, "TrimNonNumeric", true);
    A{ii} = readmatrix([pn_csv,fn_inp{ii}],opts);
end

clear opts

% Metadata & formatting
cbar_label = {'$X_{Fe}\;[at\%]$','$X_{P}\;[at\%]$','$X_{Fe}/X_{P}$','$X_{Fe}\;[at\%]$','$X_{P}\;[at\%]$','$X_{Fe}/X_{P}$'};
caxis_lims = [0,40; 0, 40; 1, 2.5;0,12; 0, 12; 1, 2.5];

%% plot

figure;
cmap = colormap(jet);
cmap(1,:) = [0,0,0]; % set first color to black
colormap(cmap);

for ii=1:3
    ax{ii} = subplot(2,3,ii);
   
    imagesc(A{ii});
    cbar{ii} = colorbar;
    cbar{ii}.Label.Interpreter = 'latex';
    cbar{ii}.TickDirection = 'out';
    cbar{ii}.Label.FontSize = 14;
    cbar{ii}.Label.String = cbar_label{ii};
    ax{ii}.XTick = [];
    ax{ii}.YTick = [];
    axis equal tight
    caxis(ax{ii},caxis_lims(ii,:));
    
    ax{ii+3} = subplot(2,3,ii+3);
    imagesc(A{ii+3});
    cbar{ii+3} = colorbar;
    cbar{ii+3}.Label.Interpreter = 'latex';
    cbar{ii+3}.TickDirection = 'out';
    cbar{ii+3}.Label.FontSize = 14;
    cbar{ii+3}.Label.String = cbar_label{ii};
    ax{ii+3}.XTick = [];
    ax{ii+3}.YTick = [];
    axis equal tight
    caxis(ax{ii+3},caxis_lims(ii+3,:));
end

if save_figure
    saveas(gcf,[pn_out,fn_out],'epsc');
    disp(['Saved to file:',pn_out,fn_out] );
else
    disp('Figure not saved!');
end
