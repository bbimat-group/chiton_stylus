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
% This script plots data underlying Figure 5
%
% This script was written by: Maya Kompella and Derk Joester
%
% This script makes use of functions distributed as part of AJ Johnson's error_ellipse.m
% AJ Johnson (2021). error_ellipse
% (https://www.mathworks.com/matlabcentral/fileexchange/4705-error_ellipse),
% MATLAB Central File Exchange. Retrieved February 27, 2021.


%% clean up
close all
clear variables
clc

% Import data from text files
fn_inp = {'Fig5_Core_dry_processed.csv','Fig5_Stylus_longitudinal_dry.csv','Fig5_Stylus_transversal_wet.csv'};
fn_out = {'Fig 5G'};

%% set current folder to folder containing the script
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

% give paths to data and for output relative to current folder
pn_csv = '../csv/';
pn_out = '../';

% flags
save_figure = false; 

% confidence interval for ellipse
conf = 0.95;

% colors
cstr = ["r","g","b"];

% legends
expt = ["core, dry","stylus, dry","stylus, wet"];

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
varnames = ["E_r","H"];
opts1.VariableNames = varnames;
opts1.VariableTypes = ["double", "double"];
opts2.VariableNames = varnames;
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

%% determine confidence ellipses

% calculate covariance matrix and center of distribution on log axes
% code based on AJ Johnson's error_ellipse (see below)
for ii = 1:length(fn_inp)
    lgE = log10(T{ii}.E_r);
    lgH = log10(T{ii}.H);
    mu_lgE  = mean(lgE);
    mu_lgH  = mean(lgH);
    C   = cov(lgE,lgH);

    [r,c] = size(C);
    k = sqrt(qchisq(conf,r));

    n=100; % Number of points around ellipse
    chi=0:pi/n:2*pi; % angles around a circle

    [eigvec,eigval] = eig(C); % Compute eigen-stuff
    xy(:,:,ii) = 10.^([mu_lgE,mu_lgH] + k*[cos(chi'),sin(chi')] * sqrt(eigval) * eigvec'); % Transformation
end


%% plot data for ashby plot in figure 5
figure('Position',[0 0 1000 500]);
for ii=1:length(fn_inp)
    % loglog(T{ii}.E_r,T{ii}.H,cstr(ii)+".",'MarkerSize',14,'DisplayName',expt(ii));hold on;
    % patch(xy(:,1,ii),xy(:,2,ii),cstr(ii),'FaceAlpha',0.2,'EdgeColor','None','DisplayName',[num2str(conf*100,'%2.0f'),'% conf. ellipse']);
    loglog(xy(:,1,ii),xy(:,2,ii),cstr(ii),'DisplayName',[num2str(conf*100,'%2.0f'),'% conf. ellipse for ',char(expt(ii))]);hold on;
end
%axis square;
ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k';
ax1.TickDir = 'out';
ax1.XLim = [1,100];
ax1.YLim = [0.1,10];
ylabel([T{1}.Properties.VariableNames{2},' [',T{1}.Properties.VariableUnits{2},']'])
xlabel([T{1}.Properties.VariableNames{1},' [',T{1}.Properties.VariableUnits{1},']'])
legend('Location','best')


%%
if save_figure
    saveas(gcf,[pn_out,fn_out{1}],'epsc');
    disp(['Saved to file:',pn_out,fn_out{1}] );
else
    disp('Figure not saved!');
end

%% the functions below are from AJ Johnson's error_ellipse.m
% AJ Johnson (2021). error_ellipse
% (https://www.mathworks.com/matlabcentral/fileexchange/4705-error_ellipse),
% MATLAB Central File Exchange. Retrieved February 27, 2021.

function x=qchisq(P,n)
% QCHISQ(P,N) - quantile of the chi-square distribution.
if nargin<2
  n=1;
end

s0 = P==0;
s1 = P==1;
s = P>0 & P<1;
x = 0.5*ones(size(P));
x(s0) = -inf;
x(s1) = inf;
x(~(s0|s1|s))=nan;

for ii=1:14
  dx = -(pchisq(x(s),n)-P(s))./dchisq(x(s),n);
  x(s) = x(s)+dx;
  if all(abs(dx) < 1e-6)
    break;
  end
end
end
%---------------------------------------------------------------
function F=pchisq(x,n)
% PCHISQ(X,N) - Probability function of the chi-square distribution.
if nargin<2
  n=1;
end
F=zeros(size(x));

if rem(n,2) == 0
  s = x>0;
  k = 0;
  for jj = 0:n/2-1;
    k = k + (x(s)/2).^jj/factorial(jj);
  end
  F(s) = 1-exp(-x(s)/2).*k;
else
  for ii=1:numel(x)
    if x(ii) > 0
      F(ii) = quadl(@dchisq,0,x(ii),1e-6,0,n);
    else
      F(ii) = 0;
    end
  end
end
end
%---------------------------------------------------------------
function f=dchisq(x,n)
% DCHISQ(X,N) - Density function of the chi-square distribution.
if nargin<2
  n=1;
end
f=zeros(size(x));
s = x>=0;
f(s) = x(s).^(n/2-1).*exp(-x(s)/2)./(2^(n/2)*gamma(n/2));
end
