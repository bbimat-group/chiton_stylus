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
% This script plots data underlying Figure S16
%
% This script was written by: Maya Kompella and Derk Joester

%% Note
% the latest version of this figure was plotted using a python script
%% Clean up
close all
clear variables
clc

%% Declare constants
fn_inp = {'FigS16_santabarabaraite_processed.csv','FigS16_stylus_processed.csv'};
fn_out = 'FigS16';

% Offset for plots; first plot on bottom
offset = [0,0.3];

% linestyle for plots
f={'k','r'};

% legend entries
legendstr = {'SBB','Stylus'};
% Flags for output
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
opts1 = delimitedTextImportOptions("NumVariables", 6);  % for data import
opts2 = opts1; % for unit import

% Specify range and delimiter
opts1.DataLines = [3, Inf];
opts1.Delimiter = ",";
opts2.DataLines = [2, 2];
opts2.Delimiter = opts1.Delimiter;

% Specify column names and types
opts1.VariableNames = ["Energy", "NormalizedIntensity"];
opts1.VariableTypes = ["double", "double"];
opts2.VariableNames = opts1.VariableNames;
opts2.VariableTypes = ["string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
for ii=1:length(fn_inp)
    T{ii} = readtable([pn_csv,fn_inp{ii}],opts1);
    U     = readtable([pn_csv,fn_inp{ii}],opts2);
    T{ii}.Properties.VariableUnits = U{:,:};
end

%% Plot data
for ii = 2:-1:1
   %plot and shift by offset
   plot(T{ii}.Energy,T{ii}.NormalizedIntensity+offset(ii),'LineWidth',1,'color',f{ii},'DisplayName',legendstr{ii}); hold on; 

    ax1 = gca;
    ax1.XLim = [7050 7250];
    ax1.TickDir = 'out';
    ax1.YTick = [];
    ylabel(['Normalized intensity [',T{ii}.Properties.VariableUnits{2},']'],'color','black');
    xlabel(['Energy [',T{ii}.Properties.VariableUnits{1},']']); 
    ax1.XGrid = 'on';
    ax1.XMinorTick = 'on';
    ax1.XMinorGrid = 'on';

    ax1.GridAlpha = 0.8;
    ax1.MinorGridAlpha = 0.8;
    legend('location','best')

end


if save_figure==true
    saveas(gcf,[pn_out,fn_out],'epsc');
    disp(['Saved to file:',pn_out,fn_out] );
else
    disp('Figure not saved!');
end   
