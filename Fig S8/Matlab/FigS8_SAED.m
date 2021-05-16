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
% This script plots data underlying Panels W, X, and Y of Figure S8
%
% This script was written by: Maya Kompella and Derk Joester

%% clean up
clear variables
close all
clc

%% Declare constants
fn = {'SAED_Radial average_mineralized_Stylus.csv','Fig_S8W.tif'};
fn_out = 'Fig S8W2Y';

% flags for output
save_figure = false; 

%% set current folder to folder containing the script
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

% give paths to data and for output relative to current folder
pn_csv = '../csv/';
pn_out = '../';


% % read image
I = imread([pn_out,fn{2}]);
[Nrow,Ncol] = size(I);

% import radial integral generated in GSM
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [3, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["frequency", "Intensity"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
RInt = readtable([pn_csv,fn{1}], opts);

% scale for SAED pattern
scale_old = 0.004; % [nm^-1/pixel]

% beamcenter
xci = 2068;
yci = 1998;

fvx = ([1:Ncol]-xci)*scale_old;
fvy = ([1:Ncol]-yci)*scale_old;

%% plot SAED pattern and superpose radial integral
close all
Is = (I-mean(I,'all'))/std(I,[],'all'); 
cmax = quantile(Is(:),0.98);
cmin = quantile(Is(:),0.05);

figure('Units','inches','Position',[592 1164 1369 625]/72);
tiledlayout(2,4);
nexttile([2,2]);
imagesc(fvx,fvy,Is); hold on;
xlabel('frequency [nm^{-1}]');
ylabel('frequency [nm^{-1}]');
ind = RInt.frequency >= 0.5;
plot(RInt.frequency(ind),16*RInt.Intensity(ind)/max(RInt.Intensity(ind)),'y.');
xline(0,'y--','LineWidth',0.5);
yline(0,'y--','LineWidth',0.5);
ax = gca;
ax.TickDir = 'out';
axis xy
daspect([1,1,1]);
colormap gray
caxis([cmin,cmax]);

%% plot 
% ax1 = nexttile;
% plot(RInt.frequency,RInt.Intensity,'.');hold on;
% 
% ax1.XLim = [0.5,4];
% ax1.YLim = [0,100];
% fit in frequency domain

idx1 =  RInt.frequency>=1 & RInt.frequency<=4;
idx2 =  RInt.frequency>=1.5 & RInt.frequency<=2.5;
f_exp = RInt.frequency(idx1 & ~idx2);
I_exp = RInt.Intensity(idx1 & ~idx2);

fitfun = @(A,x) 1./(A(1)*x.^2 + A(2)*x + A(3));

[A1,rn] = lsqcurvefit(fitfun, [1 1 0], f_exp,I_exp);
f_ax = linspace(1,4,100);

I_res = RInt.Intensity-fitfun(A1,RInt.frequency);
I_res_smoothed = movmean(I_res,20);


[max_I,max_idx] = max(I_res_smoothed(idx2));
max_f = RInt.frequency(idx2);

max_idx = find(RInt.frequency == max_f(max_idx) & RInt.frequency >=2);
max_f   = RInt.frequency(max_idx);

% plot
ax2  = nexttile;
plot(RInt.frequency,RInt.Intensity,'.'); hold on;
plot(f_ax,fitfun(A1,f_ax),'r--')
ax2.XLim = [1,4];
ax2.YLim = [0,50];
ax2.TickDir = 'out';
xlabel('Frequency [nm^{-1}]')
ylabel('Intensity')

ax3 = nexttile;
plot(RInt.frequency,I_res,'g.'); hold on
plot(RInt.frequency,I_res_smoothed,'m-','LineWidth',1)
plot([max_f,max_f],[max_I,max_I+1],'r','LineWidth',1);
text(max_f,max_I+1.5,num2str(max_f,'%1.2f'),'HorizontalAlignment','center')

ax3.XLim = [1,4];
ax3.YLim = [-1,10];
xlabel('Frequency [nm^{-1}]')
ylabel('Intensity-Bkg [a.u.]')
ax3.TickDir = 'out';

%% fit in real space
RInt.Distance = 1./RInt.frequency;
idx1 =  RInt.Distance>=0.25 & RInt.Distance<=1;
idx2 =  RInt.Distance>=0.4 & RInt.Distance<=0.67;

d_exp = RInt.Distance(idx1 & ~idx2);
I_exp = RInt.Intensity(idx1 & ~idx2);

fitfun = @(A,x) A(1)*x.^2 + A(2)*x + A(3);

[A2,rn2] = lsqcurvefit(fitfun, [1 1 0], d_exp,I_exp);
d_ax = linspace(0.1,1,100);

I_res = RInt.Intensity-fitfun(A2,RInt.Distance);
I_res_smoothed = movmean(I_res,40);

[max_I,max_idx] = max(I_res_smoothed(idx2));
max_d = RInt.Distance(idx2);

max_idx = find(RInt.Distance == max_d(max_idx) & RInt.Distance <=1);
max_d   = RInt.Distance(max_idx);

% plot
ax4  = nexttile;
plot(10*d_exp,I_exp)
plot(10*RInt.Distance,RInt.Intensity,'.'); hold on;
plot(d_ax*10,fitfun(A1,d_ax),'r--')
ax4.XLim = [2.5,7];
ax4.TickDir = 'out';
%ax1.YLim = [0,50];
xlabel('Distance [\AA]','interpreter','latex')
ylabel('Intensity')
 
ax5 = nexttile;
plot(RInt.Distance*10,RInt.Intensity-fitfun(A1,RInt.Distance),'g.'); hold on;
plot(RInt.Distance*10,I_res_smoothed,'m-','LineWidth',1)
plot(10*[max_d,max_d],[max_I,max_I+1],'r','LineWidth',1);
text(10*max_d,max_I+1.5,num2str(10*max_d,'%1.2f'),'HorizontalAlignment','center')
ax5.TickDir = 'out';

ax5.XLim = [2.5,7];
ax5.YLim = [-1,10];
xlabel('Distance [\AA]','interpreter','latex')
ylabel('Intensity-Bkg [a.u.]')

if save_figure
    saveas(gcf,[pn_out,fn_out],'epsc');
    disp(['Saved to file:',pn_out,fn_out] );
else
    disp('Figure not saved!');
end