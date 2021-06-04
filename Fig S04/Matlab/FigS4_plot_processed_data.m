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
% This script plots data underlying Figure S4
%
% This script was written by: Maya Kompella and Derk Joester

%% clean up
close all
clear variables
clc


%% Declare constants

expt   = {'Stylus','Stylus demin.','AFP', 'ChitosanOAc','ChitosanOAc/AFP 75wt%'}
fn_inp = {'FigS4_processed_AFP.csv','FigS4_processed_Stylus_demin.csv','FigS4_processed_Stylus.csv','FigS4_processed_ChitosanOAc.csv','FigS4_processed_Composite75.csv'};
fn_out = 'FigS4';

% offset for plots here; first plot on top
offset = [4,3,2,1,0];

% linestyle for plots
f={'b','g','r','#800080','k'};

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
% Setup the Import Options and import the data
opts1 = delimitedTextImportOptions("NumVariables", 2);  % for data import
opts2 = opts1; % for unit import

% Specify range and delimiter
opts1.DataLines = [3, Inf];
opts1.Delimiter = ",";
opts2.DataLines = [2, 2];
opts2.Delimiter = opts1.Delimiter;

% Specify column names and types
opts1.VariableNames = ["Wavenumber", "Transmission"];
opts1.VariableTypes = ["double", "double"];
opts2.VariableNames = opts1.VariableNames;
opts2.VariableTypes = ["string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

for ii=1:length(fn_inp)
    T{ii} = readtable([pn_csv,fn_inp{ii}],opts1);
    U     = readtable([pn_csv,fn_inp{ii}],opts2);
    T{ii}.Properties.VariableUnits = U{:,:};
end

clear opts1 opts2 U


%% Plot data
for ii = 1:5
   %plot and shift by offset
    plot(T{ii}.Wavenumber,T{ii}.Transmission+offset(ii),'LineWidth',1,'color',f{ii}); hold on; 
end

ax1 = gca;
ax1.XLim = [600 4000];
ax1.TickDir = 'out';
ax1.YTick = [];
ylabel(['Transmission [',T{ii}.Properties.VariableUnits{2},']'],'color','black');
xlabel(['Wavenumber [',T{ii}.Properties.VariableUnits{1},']']); 
ax1.XGrid = 'on';
ax1.XMinorTick = 'on';
ax1.XMinorGrid = 'on';
set(gca, 'XDir','reverse')
ylim([-1 5.5])

ax1.GridAlpha = 0.8;
ax1.MinorGridAlpha = 0.8;
legend(expt,'location','best')

if save_figure==true
    saveas(gcf,[pn_out,fn_out],'epsc');
    disp(['Saved to file:',pn_out,fn_out] );
else
    disp('Figure not saved!');
end
