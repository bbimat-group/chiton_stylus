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
% This script plots data underlying Figure S17
%
% This script was written by: Maya Kompella and Derk Joester

%% clean up
close all
clear variables
clc


%% Declare constants

fn_inp = 'FigS17_processed.csv';
fn_out = 'FigS17';

%Set offset and color for plots
c      = {'.k'}; % << set color for plots here

% Flags for output
save_figure = false; % << keep this false unless you are sure you want to overwrite the current figure

%% set current folder to folder containing the script
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

% give paths to data and for output relative to current folder
pn_csv = '../csv/';
pn_out = '../';


%% Import data from text files
% Setup the Import Options and import the data
opts1 = delimitedTextImportOptions("NumVariables", 3);  % for data import
opts2 = opts1; % for unit import

% Specify range and delimiter
opts1.DataLines = [3, Inf];
opts1.Delimiter = ",";
opts2.DataLines = [2, 2];
opts2.Delimiter = opts1.Delimiter;

% Specify column names and types
opts1.VariableNames = ["Depth", "Load"];
opts1.VariableTypes = ["double", "double"];
opts2.VariableNames = opts1.VariableNames;
opts2.VariableTypes = ["string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";


T = readtable([pn_csv,fn_inp],opts1);
U = readtable([pn_csv,fn_inp],opts2);
T.Properties.VariableUnits = U{:,:};

clear opts1 opts2 U

%% Plot data

plot(T.Depth,T.Load,c{1},'LineWidth',1);

ax1 = gca;
ax1.XLim = [0 400];
ax1.YLim = [min(T.Load) 1500];
ylabel(['load [',T.Properties.VariableUnits{2},'])'],'color','black');
xlabel(['depth [',T.Properties.VariableUnits{1},']']);
ax1.TickDir = 'out';
ax1.XGrid = 'on';

if save_figure
    saveas(gcf,[pn_out,fn_out],'epsc');
    disp(['Saved to file:',pn_out,fn_out] );
else
    disp('Figure not saved!');
end
