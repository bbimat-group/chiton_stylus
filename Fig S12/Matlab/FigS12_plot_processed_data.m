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
% This script plots data underlying Figure S12
%
% This script was written by: Maya Kompella and Derk Joester

%% clean up
close all
clear variables
clc


%% Declare constants

fn_inp = {'FigS12_Stylus_st1.csv','FigS12_Stylus_st2.csv','FigS12_Stylus_st3.csv'};
fn_out = 'Fig S12_v2';

% offset for plots here; first plot on top
offset = [1,0.5,0]; 

% plot window
tmin = 25; %[ns]
tmax = 120;%[ns]

% linestyle for plots
f={'r.','b.','k.'};

legendnames={'st1','st2 ','st3'};

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
opts1.VariableNames = ["Time", "Intensity"];
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
close all
figure;
for ii = 1:3
   % normalize range in window to 1 decade
   idx = (T{ii}.Time >= tmin) & (T{ii}.Time<=tmax) & T{ii}.Intensity>=5;
   I = T{ii}.Intensity(idx);
   t = T{ii}.Time(idx);

   semilogy(t,I,f{ii},'LineWidth',1,'DisplayName',legendnames{ii},'MarkerSize',14);

   hold on

   ax1 = gca;
   ax1.XLim = [tmin tmax];
   ax1.YLim = [5,2e3];
   ax1.TickDir = 'out';
   xlabel(['Time [', T{ii}.Properties.VariableUnits{1}, ']']); 
   ylabel(['Intensity [' , T{ii}.Properties.VariableUnits{2}, ']'])    
   legend('location','best');
end   

if save_figure
    saveas(gcf,[pn_out,fn_out],'epsc');
    disp(['Saved to file:',pn_out,fn_out] );
else
    disp('Figure not saved!');
end
