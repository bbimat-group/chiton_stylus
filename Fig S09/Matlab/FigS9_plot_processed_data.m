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
% This script plots data underlying Figure S9
%
% This script was written by: Maya Kompella and Derk Joester

%% clean up
clear variables
close all
clc

%% Declare constants
fn = {'Fig S9C.tif','Radial_Integral_SAED.csv','Radial_Integral_SAED_Fiji.csv'};
fn_out = 'Fig S9C2G';

% flags for output
save_figure = false; 

%% set current folder to folder containing the script
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

% give paths to data and for output relative to current folder
pn_csv = '../csv/';
pn_out = '../';


% read image
I2 = imread([pn_out,fn{1}]);

% import radial integral generated in GSM
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [4, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Frequency", "Distance", "Intensity"];
opts.VariableTypes = ["double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
RInt_I2 = readtable([pn_csv,fn{2}], opts);

% import radial integral generated in Fiji
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Pos", "Intensity"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
RInt_I2_Fiji = readtable([pn_csv,fn{3}], opts);

% Clear temporary variables
clear opts

%% plot SAED

ori = [275,334];
close all
figure('Position',[200 897 1575 693]);
ax1 = subplot(2,4,[1,2,5,6]);
cmap = colormap(gray);
image(I2);hold on;
profile = ori(2)-800*rescale(RInt_I2_Fiji.Intensity);
profile(profile<0) = 0;
plot(RInt_I2_Fiji.Pos+ori(1),profile,'y','LineWidth',1);
axis off equal tight

%% fit in frequency domain

idx1 =  RInt_I2.Frequency>=1 & RInt_I2.Frequency<=5;
idx2 =  RInt_I2.Frequency>=2 & RInt_I2.Frequency<=4;
f_exp = RInt_I2.Frequency(idx1 & ~idx2);
I_exp = RInt_I2.Intensity(idx1 & ~idx2);

fitfun = @(A,x) 1./(A(1)*x.^2 + A(2)*x + A(3));

[A1,rn] = lsqcurvefit(fitfun, [1 1 0], f_exp,I_exp);
f_ax = linspace(1,5,100);

I_res = RInt_I2.Intensity-fitfun(A1,RInt_I2.Frequency);
I_res_smoothed = movmean(I_res,20);


[max_I,max_idx] = max(I_res_smoothed(idx2));
max_f = RInt_I2.Frequency(idx2);

max_idx = find(RInt_I2.Frequency == max_f(max_idx) & RInt_I2.Frequency >=2);
max_f   = RInt_I2.Frequency(max_idx);

% plot
ax2  = subplot(2,4,3);
plot(RInt_I2.Frequency,RInt_I2.Intensity,'.'); hold on;
plot(f_ax,fitfun(A1,f_ax),'r--')
ax2.XLim = [1.2,5];
ax2.YLim = [0,50];
ax2.TickDir = 'out';
xlabel('Frequency [nm^{-1}]')
ylabel('Intensity')

ax3 = subplot(2,4,7);
plot(RInt_I2.Frequency,I_res,'g.'); hold on
plot(RInt_I2.Frequency,I_res_smoothed,'m-','LineWidth',1)
plot([max_f,max_f],[max_I,max_I+1],'r','LineWidth',1);
text(max_f,max_I+1.5,num2str(max_f,'%1.2f'),'HorizontalAlignment','center')

ax3.XLim = [1.2,5];
ax3.YLim = [-1,10];
xlabel('Frequency [nm^{-1}]')
ylabel('Intensity-Bkg [a.u.]')
ax3.TickDir = 'out';

%% fit in real space
idx1 =  RInt_I2.Distance>=0.2 & RInt_I2.Distance<=0.9;
idx2 =  RInt_I2.Distance>=0.25 & RInt_I2.Distance<=0.5;

d_exp = RInt_I2.Distance(idx1 & ~idx2);
I_exp = RInt_I2.Intensity(idx1 & ~idx2);

fitfun = @(A,x) A(1)*x.^2 + A(2)*x + A(3);

[A1,rn] = lsqcurvefit(fitfun, [1 1 0], d_exp,I_exp);
d_ax = linspace(0.1,1,100);

I_res = RInt_I2.Intensity-fitfun(A1,RInt_I2.Distance);
I_res_smoothed = movmean(I_res,20);

[max_I,max_idx] = max(I_res_smoothed(idx2));
max_d = RInt_I2.Distance(idx2);

max_idx = find(RInt_I2.Distance == max_d(max_idx) & RInt_I2.Distance <=1);
max_d   = RInt_I2.Distance(max_idx);

% plot
ax4  = subplot(2,4,4);
plot(10*d_exp,I_exp)
plot(10*RInt_I2.Distance,RInt_I2.Intensity,'.'); hold on;
plot(d_ax*10,fitfun(A1,d_ax),'r--')
ax4.XLim = [2,7];
ax4.TickDir = 'out';
%ax1.YLim = [0,50];
xlabel('Distance [\AA]','interpreter','latex')
ylabel('Intensity')
 
ax5 = subplot(2,4,8);
plot(RInt_I2.Distance*10,RInt_I2.Intensity-fitfun(A1,RInt_I2.Distance),'g.'); hold on;
plot(RInt_I2.Distance*10,I_res_smoothed,'m-','LineWidth',1)
plot(10*[max_d,max_d],[max_I,max_I+1],'r','LineWidth',1);
text(10*max_d,max_I+1.5,num2str(10*max_d,'%1.2f'),'HorizontalAlignment','center')
ax5.TickDir = 'out';

ax5.XLim = [2,7];
ax5.YLim = [-1,10];
xlabel('Distance [\AA]','interpreter','latex')
ylabel('Intensity-Bkg [a.u.]')

if save_figure
    saveas(gcf,[pn_out,fn_out],'epsc');
    disp(['Saved to file:',pn_out,fn_out] );
else
    disp('Figure not saved!');
end