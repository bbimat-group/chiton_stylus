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
% This script plots data underlying Figure S14
%
% This script was written by: Maya Kompella and Derk Joester

%% clean up
close all
clear variables
clc

%% Declare constants

% Import data from text files
fn_inp = {'FigS14a_processed_core_at_54K.csv',...
          'FigS14a_processed_core_at_54K_fit.csv',...
          'FigS14a_processed_core_at_RT.csv',...
          'FigS14a_processed_core_at_RT_fit.csv',...
          'FigS14c_processed_st1_at_54K.csv',...
          'FigS14c_processed_st1_at_54K_fit.csv',...
          'FigS14c_processed_st1_at_RT.csv',...
          'FigS14c_processed_st1_at_RT_fit.csv',...
          'FigS14b_processed_core_at_54K.csv',...
          'FigS14b_processed_core_at_RT.csv',...
          'FigS14d_processed_st1_at_54K.csv',...
          'FigS14d_processed_st1_at_RT.csv'};
      
fn_out = 'Fig S14_v2';

% offset for plots here; first plot on top
offset = 0.2+[0,0,1,1,0,0,0.5,0.5,0,.75,0,.75];

% linestyle for plots
c={'r','r','k','k','r','r','k','k','r','k','r','k'};

%Legend names
legendnames={
          'core @ 54K',...
          'core @ 54K_fit',...
          'core @ RT',...
          'core @ RT_fit',...
          'st1 @ 54K',...
          'st1 @ 54K_fit',...
          'st1 @ RT',...
          'st1 @ RT_fit',...
          'core @ 54K',...
          'core @ RT',...
          'st1 @ 54K',...
          'st1 @ RT'};

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
opts1 = delimitedTextImportOptions("NumVariables", 2);  % for data import
opts2 = opts1; % for unit import

% Specify range and delimiter
opts1.DataLines = [3, Inf];
opts1.Delimiter = ",";
opts2.DataLines = [2, 2];
opts2.Delimiter = opts1.Delimiter;

% Specify column names and types
opts1.VariableNames = ["Time", "NormLogIntensity"];
opts1.VariableTypes = ["double", "double"];
opts2.VariableNames = opts1.VariableNames;
opts2.VariableTypes = ["string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

for ii=1:8
    T{ii} = readtable([pn_csv,fn_inp{ii}],opts1);
    U     = readtable([pn_csv,fn_inp{ii}],opts2);
    T{ii}.Properties.VariableUnits = U{:,:};
end

% Setup the Import Options and import the data
opts1 = delimitedTextImportOptions("NumVariables", 2);  % for data import
opts2 = opts1; % for unit import

% Specify range and delimiter
opts1.DataLines = [3, Inf];
opts1.Delimiter = ",";
opts2.DataLines = [2, 2];
opts2.Delimiter = opts1.Delimiter;

% Specify column names and types
opts1.VariableNames = ["Velocity","NormIntensity"];
opts1.VariableTypes = ["double", "double"];
opts2.VariableNames = opts1.VariableNames;
opts2.VariableTypes = ["string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

for ii=9:12
    T{ii} = readtable([pn_csv,fn_inp{ii}],opts1);
    U     = readtable([pn_csv,fn_inp{ii}],opts2);
    T{ii}.Properties.VariableUnits = U{:,:};
end

clear opts1 opts2 U


%% Plot
figure('Units','inches','Position',[8 8 14.5 8+6.5/4*3])
subplot(2,2,1)
for ii = 1:2:3
    
    % interpolate fit on data time axis and calculate residual
    residual = 10.^T{ii}.NormLogIntensity-interp1(T{ii+1}.Time,10.^T{ii+1}.NormLogIntensity,T{ii}.Time);
    % plot
    plot(T{ii}.Time,T{ii}.NormLogIntensity+offset(ii),'.','color',c{ii},'HandleVisibility','off'); hold on
    plot(T{ii+1}.Time,T{ii+1}.NormLogIntensity+offset(ii),'LineWidth',1,'color',c{ii},'DisplayName',legendnames{ii}); hold on; 
    lgd = legend();
    lgd.NumColumns = 1;
    ax1 = gca;
    ax1.XLim = ([25 120]);
    %ax1.YLim = ([0 2]);
    set(ax1, 'box','on','YTickLabel',[],'YTick',[])
    ylabel('normalized log intensity [a.u.]','color','black');
    xlabel(['time [',T{ii}.Properties.VariableUnits{1},']']); 
    ax1.TickDir = 'out';
    grid on
end

subplot(2,2,3)
for ii = 5:2:7
    plot(T{ii}.Time,T{ii}.NormLogIntensity+offset(ii),'.','color',c{ii-4},'HandleVisibility','off'); hold on
    plot(T{ii+1}.Time,T{ii+1}.NormLogIntensity+offset(ii),'LineWidth',1,'color',c{ii},'DisplayName',legendnames{ii}); hold on; 
    lgd = legend;
    lgd.NumColumns = 1;
    ax1 = gca;
    ax1.XLim = ([25 120]);
    ax1.YLim = ([0 1.5]);
    set(ax1, 'box','on','YTickLabel',[],'YTick',[])
    ylabel(['normalized log intensity [',T{ii}.Properties.VariableUnits{2},']'],'color','black');
    xlabel(['time [',T{ii}.Properties.VariableUnits{1},']']); 
    ax1.TickDir = 'out';
    grid on
end

subplot(2,2,2)
for ii=9:10
    plot(T{ii}.Velocity,T{ii}.NormIntensity+offset(ii),'LineWidth',1,'color',c{ii},'DisplayName',legendnames{ii}); hold on
    lgd = legend;
    lgd.NumColumns = 1;
    ax1 = gca;
    ax1.XLim = ([-5 5]);
    ax1.YLim = ([0 2.5]);
    set(ax1, 'box','on','YTickLabel',[],'YTick',[])
    ylabel(['normalized intensity [',T{ii}.Properties.VariableUnits{2},']'],'color','black');
    xlabel(['velocity [',T{ii}.Properties.VariableUnits{1},']']); 
    ax1.TickDir = 'out';
    grid on
end

subplot(2,2,4)
for ii=11:12
    plot(T{ii}.Velocity,T{ii}.NormIntensity+offset(ii),'LineWidth',1,'color',c{ii},'DisplayName',legendnames{ii}); hold on
    lgd = legend;
    lgd.NumColumns = 1;
    ax1 = gca;
    ax1.XLim = ([-5 5]);
    ax1.YLim = ([0 2.5]);
    set(ax1, 'box','on','YTickLabel',[],'YTick',[])
    ylabel('normalized intensity [a.u.]','color','black');
    xlabel(['velocity [',T{ii}.Properties.VariableUnits{1},']']); 
    ax1.TickDir = 'out';
    grid on
end

if save_figure==true
    saveas(gcf,[pn_out,fn_out],'epsc');
    disp(['Saved to file:',pn_out,fn_out] );
else
    disp('Figure not saved!');
end   