%% 
% This script is part of a collection used to visualize data associated
% with a manuscript published in PNAS entitled:
%
% "Persistent polyamorphism in the chiton tooth: 
% from a new biomineral to inks for additive manufacturing."
% 
% by Linus Stegbauer, E. Ercan Alp, Paul J. M. Smeets, Robert Free, Shay G. Wallace, Mark C. Hersam, and Derk Joester* 
%
% DOI: https://doi.org/10.1073/pnas.2020160118
% 
% This script plots data underlying panels E and F of Figure S6
%
% This script was written by: Maya Kompella and Derk Joester

%% clean up
close all
clear variables
clc


%% Declare constants

fn_inp = {'FigS6_e_processed.csv','FigS6_f_processed.csv'};
fn_out = 'Fig S6EF';

% flags for output
save_figure = false; 

%% set current folder to folder containing the script
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

% give paths to data and for output relative to current folder
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
cbar_label = {'Mineral fraction wt%','Mineral fraction vol%'};
caxis_lims = [0,30; 0, 30];

%% plot

figure;
cmap = colormap(jet);
cmap(1,:) = [0,0,0]; % set first color to black
cmap(2,:) = [0,0,0];
colormap(cmap);

for ii=1:2
    ax{ii} = subplot(1,2,ii);
   
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
    
end

if save_figure
    saveas(gcf,[pn_out,fn_out],'epsc');
    disp(['Saved to file:',pn_out,fn_out] );
else
    disp('Figure not saved!');
end
