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
% This script plots data underlying Panels B and C of Figure S18
%
% This script was written by: Maya Kompella and Derk Joester

%% clean up
close all
clear variables
clc

% Import data from text files
fn_inp = {'FigS18_transverse_processed.csv','FigS18_trailing_edge_processed.csv'};
fn_out = {'Fig S18BC'};

%% set current folder to folder containing the script
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

% give paths to data and for output relative to current folder
pn_csv = '../csv/';
pn_out = '../';

% flags
save_figure = true; 

%% Import data from text files
% Setup the Import Options and import the data
opts1 = delimitedTextImportOptions("NumVariables", 5);  % for data import
opts2 = opts1; % for unit import

% Specify range and delimiter
opts1.DataLines = [2, Inf];
opts1.Delimiter = ",";
opts2.DataLines = [2, 2];
opts2.Delimiter = opts1.Delimiter;

% Specify column names and types
varnames = ["Distance","E_r","H","stdev_E","stdev_H"];
opts1.VariableNames = varnames;
opts1.VariableTypes = ["double", "double","double", "double","double"];
opts2.VariableNames = varnames;
opts2.VariableTypes = ["string", "string","string", "string","string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

for ii=1:length(fn_inp)
    T{ii} = readtable([pn_csv,fn_inp{ii}],opts1);
    U     = readtable([pn_csv,fn_inp{ii}],opts2);
    T{ii}.Properties.VariableUnits = U{:,:};
end

clear opts1 opts2 U


%% plot figure 17
figure('Position',[0 0 1000 500]);
subplot(1,2,1)
errorbar(T{1}.Distance,T{1}.E_r,T{1}.stdev_E,'^k','MarkerFaceColor','k')
box off
ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k';
ax1.TickDir = 'out';
ylim([5 35])
xlim([-5 250])
ylabel(['Red. Modulus [',T{1}.Properties.VariableUnits{2},']'])
xlabel(['Distance [',T{1}.Properties.VariableUnits{1},']'])


ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');
ax2.XColor = 'none';
ax2.YColor = 'r';
ax2.XTick = [];
ax2.TickDir = 'out';
ylabel(['Hardness [',T{1}.Properties.VariableUnits{3},']'],'color','r')

hold on
errorbar(T{1}.Distance,T{1}.H,T{1}.stdev_H,'sr','MarkerFaceColor','r')
ylim([0.3 2])
xlim([-5 250])
title('Linescan transverse')

%% plot figure 18
subplot(1,2,2)
errorbar(T{2}.Distance,T{2}.E_r,T{2}.stdev_E,'^k','MarkerFaceColor','k')
box off
ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k';
ax1.TickDir = 'out';
ylim([5 35])
xlim([-5 250])
ylabel(['Red. Modulus [',T{2}.Properties.VariableUnits{2},']'])
xlabel(['Distance [',T{2}.Properties.VariableUnits{1},']'])


ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');
ax2.XColor = 'none';
ax2.YColor = 'r';
ax2.XTick = [];
ax2.TickDir = 'out';
ylabel(['Hardness [',T{2}.Properties.VariableUnits{3},']'],'color','r')

hold on
errorbar(T{2}.Distance,T{2}.H,T{2}.stdev_H,'sr','MarkerFaceColor','r')
ylim([0.3 2])
xlim([-5 250])
title('Linescan trailing edge')
%%
if save_figure
    saveas(gcf,[pn_out,fn_out{1}],'epsc');
    disp(['Saved to file:',pn_out,fn_out{1}] );
else
    disp('Figure not saved!');
end