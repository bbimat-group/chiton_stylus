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
% This script plots data underlying Panels N to U of Figure S8 
%
% This script was written by: Maya Kompella and Derk Joester

%% Clean up
close all
clear variables
clc

%% Set current folder to folder containing the script
mfile_name     = mfilename('fullpath');
[pn,name,ext]  = fileparts(mfile_name);
cd(pn);

% give paths to data and for output relative to current folder
pn_csv = '../csv/';
pn_out = '../';
pn_in  = '../../Fig S7/'

%% Declare constants 

% Filenames
fn_in = 'Fig_S7A.tif';
fn_out = 'Fig_S8N2U.eps';

% flags for output
save_figure = false; 

% Parameters for TEM images
scale = 0.202; % nm/pixel
Nbig  = 4096;

% define ROIs in physical coordinates
xlims = [423+[0,25]; 681 + [0,25]; 387 + [0,25]; 607 + [0,25]]; %[nm]
ylims = [398+[0,25]; 114 + [0,25]; 348 + [0,25]; 679 + [0,25]]; %[nm]

% define point where line profiles cross
xhair = [[433.09,415.31];[693.452,128.872];[395.71,364.402];[617.904,693.864]]; %[nm]

% declare functions
butter_fun  = @(D,cutoff,order) 1./(1 + (D./cutoff).^(2*order));

%% Process image
% read image
I  = double(imread([pn_in,filesep,fn_in]));
[rows,cols] = size(I);

% pad original image to 2^12 size
Ipad = zeros(Nbig);
Ipad(4+[0:rows-1],9+[0:cols-1]) = I;

% set real space coordinates
yv = [0:Nbig-1]*scale; % [nm]
xv = [0:Nbig-1]*scale; % [nm]

% center and scale
Ipad = (Ipad-mean(Ipad,'all'))/std(Ipad,[],'all');

% DFT and magnitude
Fbig = fftshift(fft2(Ipad));
Mbig = abs(Fbig);

% generate frequency axis
if mod(Nbig,2)==0
    kbig=-Nbig/2:Nbig/2-1; % N even
else
    kbig=-(Nbig-1)/2:(Nbig-1)/2; % N odd
end

Dbig = Nbig*scale; % [nm]
fbig = kbig/Dbig;  % [nm^-1]

%% apply Butterworth lowpass filter

% generate ndgrid in frequency space
[fbY,fbX] = ndgrid(fbig,fbig);

% calculate distance from [0,0]
fDb = sqrt((fbX).^2+(fbY).^2);

% set order and cutoff for BW lowpass filter
order = 4;
lp = [1.5,1/2,1/3,1/4];

% generate filter, multiply with DFT and generate filtered image by inverse
% DFT. This takes a few seconds.

for ii = 1:length(lp)
    H_lpb = butter_fun(fDb,lp(ii),order);
    F_lpb = H_lpb.*Fbig;
    I_lpb(:,:,ii) = ifft2(ifftshift(F_lpb));
end

%% Get indices for line profiles
for ii=1:4
    % get indices from real space positions
    xind(ii) = find(xv>=xhair(ii,1),1,'first');
    yind(ii) = find(yv>=xhair(ii,2),1,'first');
    xvind{ii} = find((xv >= xlims(ii,1)) & (xv <= xlims(ii,2)));
    yvind{ii} = find((yv >= ylims(ii,1)) & (yv <= ylims(ii,2)));
end

%% plot grid of linked plots to enable exploring ROIs
% figure('Units','inches','Position',[10,10,15,15]);
% tiledlayout(4,4,'TileSpacing','compact');
% 
% ctr=0;
% clear ax
% for ii = 1:4 % interate over ROIs
%     for jj = 1:length(lp) %  iterate over cutoffs
%        ctr=ctr+1;
%        ax(ctr) = nexttile;
%        imagesc(xv,yv,squeeze(I_lpb(:,:,jj)));
%        ax(ctr).TickDir = 'out';
%        xlim(xlims(ii,:));
%        ylim(ylims(ii,:));
%        daspect([1,1,1]);
%        titlestr = "lowpass f< " + num2str(lp(jj),'%2.2f') + " order = " + num2str(order,'%1.0f');
%        title(titlestr);
%     end
%     linkaxes(ax(ctr-3:ctr),'xy');
% end
% colormap gray

%% Plot figure for supplemental
figure('Units','inches','Position',[10 3 11.6111 22.1528]);
tiledlayout(8,4,'TileSpacing','compact');
clear ax;

panelstr = ["A","B","C","D","E"];
loc = [17,19,25,27];

% plot overview - original (padded) image
ctr=1;
ax(ctr)=nexttile([4,4]);
imagesc(xv,yv,Ipad); hold on;
ax(ctr).TickDir = 'out';
daspect([1,1,1]);
colormap gray;
text(7,15,panelstr(1),'Color','w','FontSize',20);

% annotate position of ROIs in panel A
for ii = 1:4
   axes(ax(1));
   verts = combvec(xlims(ii,:),ylims(ii,:));
   patch('Vertices',verts','Faces',[1 2 4 3],'EdgeColor','w','FaceColor','none','LineWidth',1);
   text(verts(1,4)+5,verts(2,4)-6,panelstr(ii+1),'Color','w','FontSize',16); 
end

for ii = 1:4  
   % plot ROIs w/o bandpass filter
   ctr=ctr+1;
   ax(ctr) = nexttile(loc(ii));
   % give coordinates in nm relative to corner
   xv_sel = xv(xvind{ii})-xv(xvind{ii}(1));
   yv_sel = yv(yvind{ii})-yv(yvind{ii}(1));
   imagesc(xv_sel,yv_sel,Ipad(yvind{ii},xvind{ii})); 
   hold on;
   
   ax(ctr).TickDir = 'out';
   ax(ctr).XTick = [];
   daspect([1,1,1]);
   titlestr = 'unfiltered';
   title(titlestr);
   % annotate position of line profiles
   xline(xhair(ii,1)-xv(xvind{ii}(1)),'w');
   yline(xhair(ii,2)-yv(yvind{ii}(1)),'w');
   text(0.03,0.9,panelstr(ii+1),'Color','w','FontSize',20,'Units','normalized');
 
   % plot ROIs with second least stringent bandpass filter
   ctr=ctr+1;
   ax(ctr) = nexttile(loc(ii)+5);
   imagesc(xv_sel,yv_sel,squeeze(I_lpb(yvind{ii},xvind{ii},2))); hold on;
   ax(ctr).TickDir = 'out';
   ax(ctr).XTick = [];
   ax(ctr).YTick = [];
   daspect([1,1,1]);
   titlestr = "lowpass f< " + num2str(lp(2),'%2.2f') + " order = " + num2str(order,'%1.0f');
   title(titlestr);
   xline(xhair(ii,1)-xv(xvind{ii}(1)),'w');
   yline(xhair(ii,2)-yv(yvind{ii}(1)),'w');
   text(0.03,0.9,panelstr(ii+1),'Color','w','FontSize',20,'Units','normalized');
end

for ii=1:4
   % extract line profiles
   yprof_1 = Ipad(yvind{ii},xind(ii));
   xprof_1 = Ipad(yind(ii),xvind{ii});
   yprof_2 = I_lpb(yvind{ii},xind(ii),2);
   xprof_2 = I_lpb(yind(ii),xvind{ii},2);
   
   % plot line profiles in the vertical direction
   ctr=ctr+1;
   ax(ctr) = nexttile(loc(ii)+1);
   plot(yprof_1,scale*[0:length(yvind{ii})-1]'); hold on;
   plot(yprof_2,scale*[0:length(yvind{ii})-1]'); 
   yline(xhair(ii,2)-yv(yvind{ii}(1)),'k');
   ax(ctr).YAxis.Direction = 'reverse';
   ax(ctr).TickDir = 'out';
   ax(ctr).XTick = [];
   ax(ctr).YTick = [];
   ylim([0,25]);
   
   % plot ine profiles in the horizontal direction
   ctr=ctr+1;
   ax(ctr) = nexttile(loc(ii)+4);
   plot(scale*[0:length(xvind{ii})-1]',xprof_1); hold on;
   plot(scale*[0:length(xvind{ii})-1]',xprof_2); 
   xline(xhair(ii,1)-xv(xvind{ii}(1)),'k');
   xlim([0,25]);
   ax(ctr).YTick = [];
   ax(ctr).TickDir = 'out';
end

if save_figure
    saveas(gcf,[pn_out,fn_out],'epsc');
    disp(['Saved to file:',pn_out,fn_out] );
else
    disp('Figure not saved!');
end