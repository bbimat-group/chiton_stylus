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
% This script plots data underlying panel C of Figure 2
%
% This script was written by: Maya Kompella and Derk Joester

%% clean up
close all
clear variables
clc


%% Declare constants

fn_inp = {'Fig2_C_st1_exp_processed.csv','Fig2_C_core_exp_processed.csv','Fig2_C_cusp_exp_processed.csv',...
          'Fig2_C_st1_fit_processed.csv','Fig2_C_core_fit_processed.csv','Fig2_C_cusp_fit_processed.csv'};
fn_out = 'Fig2_C';

%Set offset, colors, legend entry names
offset_add = [0,.5,1];
c = {'r','k','b','r','k','b'}; 
legendnames = {'St1','Core','Cusp'};

% Flags for output
save_figure = false; % << keep this false unless you are sure you want to overwrite the current figure

% set current folder to folder containing the script
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

% give pasth to data and for output relative to current folder
pn_csv = '../csv/';
pn_out = '../';

% define time window 
tmin = 55;
tmax = 67;

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
opts1.VariableNames = ["Time", "LogIntensity"];
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


%% remove outlier and rescale cusp spectra
T{3}.LogIntensity==NaN;
T{3}.LogIntensity = rescale(T{3}.LogIntensity);

%% Plot data
figure('Units','inches','Position',[8,8,8+2,8+6])
for ii = 3:-1:1
    plot(T{ii}.Time,T{ii}.LogIntensity+offset_add(ii),'.','color',c{ii},'HandleVisibility','off','MarkerSize',14); hold on
    plot(T{ii+3}.Time,T{ii+3}.LogIntensity+offset_add(ii),'LineWidth',1,'color',c{ii},'DisplayName',legendnames{ii}); hold on; 
end
xline(tmin);
xline(tmax);

ax1 = gca;
ax1.XLim = ([25 120]);
ax1.YLim = ([-0.1 2.1]);
set(ax1, 'box','on','YTickLabel',[],'YTick',[])
ylabel(['normalized log intensity [',T{ii}.Properties.VariableUnits{2},'])'],'color','black');
xlabel(['time [',T{ii}.Properties.VariableUnits{1},']']); 
ax1.TickDir = 'out';
grid on
legend('location','southeast')
lgd = legend;
lgd.NumColumns = 1;

%% save 
if save_figure
    saveas(gcf,[pn_out,fn_out],'epsc');
    disp(['Saved to file:',pn_out,fn_out]);
else
    disp('Figure not saved!');
end
