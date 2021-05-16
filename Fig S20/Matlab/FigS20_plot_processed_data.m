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
% This script plots data underlying Figure S20
%
% This script was written by: Maya Kompella and Derk Joester

%% clean up
close all
clear variables
clc


%% Declare constants

fn_inp = {'FigS20_nanoindentation_transverse_processed.csv'};
fn_out = {'Fig S20AB','Fig S20_interpolated'};

% flags for output
save_figure = true; 

%% set current folder to folder containing the script
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

% give paths to data and for output relative to current folder
pn_csv = '../csv/';
pn_out = '../';

%% Import data from text files
% Setup the Import Options and import the data
opts1 = delimitedTextImportOptions("NumVariables", 4);  % for data import
opts2 = opts1; % for unit import

% Specify range and delimiter
opts1.DataLines = [3, Inf];
opts1.Delimiter = ",";
opts2.DataLines = [2, 2];
opts2.Delimiter = opts1.Delimiter;

% Specify column names and types
opts1.VariableNames = ["Er","H","X","Y"];
opts1.VariableTypes = ["double", "double","double", "double"];
opts2.VariableNames = opts1.VariableNames;
opts2.VariableTypes = ["string", "string","string", "string"];

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
N_indents = height(T{1});
X_pos     = unique(T{1}.X); 
Y_pos     = unique(T{1}.Y);

% set by hand 
N_rows = 15;
N_cols = 21; 

% reshape table entries into an array 
H  = reshape(T{1}.H,[N_rows,N_cols]);
Er = reshape(T{1}.Er,[N_rows,N_cols]);

% calculate range in stage coordinates
X_range = max(X_pos)-min(X_pos); % [mm] corresponding to rows!
Y_range = max(Y_pos)-min(Y_pos); % [mm] corresponding to columns!

% set up x- and y-vectors (corners) for plotting
xv = [0,X_range]*1e3; % [µm]
yv = [0,Y_range]*1e3; % [µm]

% plot w/o interpolation
f=figure;
colormap(f,parula);

ax1 = subplot(1,2,1);
imagesc(ax1,fliplr(yv),xv,H); % plot array on x and y axes. 
xlabel('x [\mum]'); 
ylabel('y [\mum]')
ax1.TickDir = 'out';
axis xy equal tight

cbar_H = colorbar('TickDirection','out');
cbar_H.Label.String = ['H [',T{ii}.Properties.VariableUnits{1},']']'; 
cbar_H.TickDirection = 'out';

ax2 = subplot(1,2,2);

imagesc(ax2,fliplr(yv),xv,Er); %plot array on x and y axes
xlabel('x [\mum]'); 
ylabel('y [\mum]')
ax2.TickDir = 'out';
axis xy equal tight

cbar_Er = colorbar('TickDirection','out');
cbar_Er.Label.String = ['E_r [',T{ii}.Properties.VariableUnits{2},']']; 
cbar_Er.TickDirection = 'out';


if save_figure==true
    saveas(gcf,[pn_out,fn_out{1}],'epsc');
    disp(['Saved to file:',pn_out,fn_out{1}] );
else
    disp('Figure not saved!');
end


% %interpolation
% 
% interpolationfactor=1; % returns 2^interpolationfactor-1 interpolated points between sample values
% 
% figure(2)
% ax1 = subplot(1,2,1);
% imagesc(ax1,fliplr(yv),xv,interp2(H,interpolationfactor)); % plot array on x and y axes. 
% xlabel('x [\mum]'); 
% ylabel('y [\mum]')
% ax1.TickDir = 'out';
% axis xy equal tight
% 
% cbar_H = colorbar('TickDirection','out');
% cbar_H.Label.String = 'H [GPa]'; 
% cbar_H.TickDirection = 'out';
% 
% ax2 = subplot(1,2,2);
% 
% imagesc(ax2,fliplr(yv),xv,interp2(Er,interpolationfactor)); %plot array on x and y axes
% xlabel('x [\mum]'); 
% ylabel('y [\mum]')
% ax2.TickDir = 'out';
% axis xy equal tight
% 
% cbar_Er = colorbar('TickDirection','out');
% cbar_Er.Label.String = 'E_r [GPa]'; 
% cbar_Er.TickDirection = 'out';
% 
% if save_figure==true
%     saveas(gcf,[pn_out,fn_out{2}],'epsc');
%     disp(['Saved to file:',pn_out,fn_out{2}] );
% else
%     disp('Figure not saved!');
% end