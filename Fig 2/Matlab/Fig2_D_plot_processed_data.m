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
% This script plots data underlying panel D of Figure 2
%
% This script was written by: Maya Kompella and Derk Joester

%% clean up
close all
clear variables
clc


%% Declare constants

fn_inp = {'Fig2_D_st1_processed.csv','Fig2_D_core_processed.csv','Fig2_D_cusp_processed.csv'};
fn_out = 'Fig2_D_revised_v2';

%Set offset and color for plots
offset        = [0,1,0]; % << offset for plots here; first plot on top
c             = {'r','k','b'}; % << set color for plots here
expt          = {'St1','Core','Cusp'};

% Flags for output
save_figure = false; % << keep this false unless you are sure you want to overwrite the current figure

% set current folder to folder containing the script
%mfile_name          = mfilename('fullpath');
%[pathstr,name,ext]  = fileparts(mfile_name);
%cd(pathstr);

% give pasth to data and for output relative to current folder
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
opts1.VariableNames = ["Velocity", "NormIntensity"];
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

% % For the original submission, this showed all spectra on the same axis
% for ii = 1:3
%    plot(T{ii}.Velocity,T{ii}.LogNormIntensity+offset(ii),'color',c{ii},'LineWidth',1); hold on; 
%     
%     ax1 = gca;
%     ax1.XLim = [-10 10];
%     ax1.TickDir ='out';
%     set(ax1, 'box','on','YTickLabel',[],'YTick',[])
%     ylabel(['log(normalized intensity [',T{ii}.Properties.VariableUnits{2},'])'],'color','black');
%     xlabel(['Velocity [',T{ii}.Properties.VariableUnits{1},']']);
%     grid on
% end
% legend({'Stylus','Core','Cusp'},'location','best')

% set plot order

% following reviewer guidance, this was split into two panels
close all
figure('Units','inches','Position',[8,8,8+1.3,8+2.14])

% plot spectrum for cusp
ax1 = subplot(5,1,[1,2]);
ii = find(matches(expt,'Cusp'));
plot(ax1,T{ii}.Velocity,T{ii}.NormIntensity+offset(ii),'color',c{ii},'LineWidth',1,'DisplayName',expt{ii}); hold on; 

ax1.XLim = [-10 10];
ax1.YLim = [-0.05,1.05];
ax1.TickDir ='out';
set(ax1, 'box','on','YTickLabel',[],'YTick',[])
ylabel(['log(normalized intensity [',T{ii}.Properties.VariableUnits{2},'])'],'color','black');
xlabel(['Velocity [',T{ii}.Properties.VariableUnits{1},']']);
grid on
legend('location','best');

% plot core and stylus (st1) spectra
ax2 = subplot(5,1,[3:5]);
ii = find(matches(expt,'Core'));
plot(ax2,T{ii}.Velocity,T{ii}.NormIntensity+offset(ii),c{ii},'LineWidth',1,'DisplayName',expt{ii});hold on;

ii = find(matches(expt,'St1'));
plot(ax2,T{ii}.Velocity,T{ii}.NormIntensity+offset(ii),c{ii},'LineWidth',1,'DisplayName',expt{ii});hold on;

xlabel(['Velocity [', T{ii}.Properties.VariableUnits{1}, ']']); 
ylabel(['normalized intensity [' , T{ii}.Properties.VariableUnits{2}, '])'])

ax2.XLim = [-1 1];
ax2.YLim = [-0.05,2.05];
%ax2.XTickLabel=[-1.0:0.2:1.0];
ax2.TickDir = 'out';
ax2.YTick   = [];
ax2.XGrid = 'on';
ax2.XMinorTick = 'on';
ax2.XMinorGrid = 'on';
%ax2.GridAlpha = 0.8;
ax2.MinorGridAlpha = 0.8;

%plot lines at minima in core spectra
ii = find(matches(expt,'Core'));
min_idx = islocalmin(T{ii}.NormIntensity);
min_val_v_core =T{ii}.Velocity(min_idx);
min_val_I_core =T{ii}.NormIntensity(min_idx);
line([min_val_v_core(1),min_val_v_core(1)],[min_val_I_core(1)-0.6,min_val_I_core(1)]+offset(ii),'LineWidth',1,'Color','k');
line([min_val_v_core(2),min_val_v_core(2)],[min_val_I_core(2)-0.6,min_val_I_core(2)]+offset(ii),'LineWidth',1,'Color','k');

ii = find(matches(expt,'St1'));
min_idx = islocalmin(T{ii}.NormIntensity);
min_val_v_stylus =T{ii}.Velocity(min_idx);
min_val_I_stylus =T{ii}.NormIntensity(min_idx);
line([min_val_v_stylus(1),min_val_v_stylus(1)],[min_val_I_stylus(1)+0.6,min_val_I_stylus(1)]+offset(ii),'LineWidth',1,'Color','r');
line([min_val_v_stylus(2),min_val_v_stylus(2)],[min_val_I_stylus(2)+0.6,min_val_I_stylus(2)]+offset(ii),'LineWidth',1,'Color','r');

legend(expt([2,1]),'location','best');

if save_figure
    saveas(gcf,[pn_out,fn_out],'epsc');
    disp(['Saved to file:',pn_out,fn_out] );
else
    disp('Figure not saved!');
end

