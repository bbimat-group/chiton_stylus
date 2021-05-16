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
% This script plots data underlying Figure S7
%
% This script was written by: Maya Kompella and Derk Joester

%% Clean up
%close all
clear variables
clc

%% Set current folder to folder containing the script
mfile_name     = mfilename('fullpath');
[pn,name,ext]  = fileparts(mfile_name);
cd(pn);

% give paths to data and for output relative to current folder
pn_out = '../';

%% Declare constants 

% Filenames
    fn_in = 'Fig_S7A.tif';
    fn_out = {'Fig_S7_p1.eps','Fig_S7_p2.eps'};
% Parameters for TEM images
   scale = 0.202; % nm/pixel

% Flags:
    save_figure = true;
% Parameters for ROI
    % upper left corner for selection for FFT
    Q_col = 1600; % pixel
    Q_row = 2300; % pixel

    % number of pixels
    N = 1024;

% Parameters for ROI2
    % upper left corner for selection for FFT
    Q_col2 = 750; % pixel
    Q_row2 = 1900; % pixel

    % number of pixels
    N2 = 2048;

%% declare Functions

% rotation matrix for counter clockwise rotation about origin by angle theta (in degrees) in 2D cartesian coords 
RotM        = @(theta) [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
gauss_model = @(A,x) A(1) + A(2).*exp(-4*log(2)*((x-A(3))/A(4)).^2);
butter_fun  = @(D,cutoff,order) 1./(1 + (D./cutoff).^(2*order));
vonMises    = @(A,x) exp(A(2)*cos(x-A(1)))./(2*pi*besseli(0,A(2)));
Lorentzian  = @(A,x) (pi*A(1).*(1+((x-A(2))/A(1)).^2)).^-1;
Gaussian    = @(A,x) (A(1)*sqrt(2*pi))^-1*exp(-0.5*(x-A(2)).^2/A(1)^2);
fitfun      = @(A,x) A(1) + A(2)*(vonMises(A(3:4),x) + vonMises([A(3)-pi,A(4)],x)) + A(5)*(Gaussian(A(6:7),x) + Gaussian([A(6),A(7)-pi],x));

%% read original TEM image and create ROI
I     = double(imread([pn_out,filesep,fn_in]));
I_sel = I(Q_row+[0:N-1],Q_col+[0:N-1]);
I_sel2 = I(Q_row2+[0:N2-1],Q_col2+[0:N2-1]);

%% Create axes in physical coords
[N_rows,N_cols] = size(I);
yv = [0:N_rows-1]*scale; % [nm]
xv = [0:N_cols-1]*scale; % [nm]

yv_sel = [0:N-1]*scale; % [nm]
xv_sel = [0:N-1]*scale; % [nm]

% generate vertices and faces for a patch to indicate selection
Verts_BoxA = scale*[[Q_col,Q_col+N,Q_col+N,Q_col];[Q_row,Q_row,Q_row+N,Q_row+N]].'; % [nm]
Faces_BoxA = [1,2,3,4];

%% FFT of selection

% FFT and shift, determine Magnitude M and Phase Angle P
F = fftshift(fft2(I_sel));
M = abs(F);
P = angle(F);

% smooth magnitude
sm_stdev    = 3;
sm_filtsize = 53;
M_s = imgaussfilt(M,sm_stdev,'FilterSize',sm_filtsize);

% generate frequency axis
if mod(N,2)==0
    k=-N/2:N/2-1; % N even
else
    k=-(N-1)/2:(N-1)/2; % N odd
end

D    = N*scale; % [nm]
freq = k/D;     % [nm^-1]

%% polar transform
fR   = linspace(0,max(freq),5000);
dchi = 0.5;
chi  = 0:dchi:360-dchi;
[fY,fX] = ndgrid(freq,freq);

W   = unwrapPolar(fX,fY,M,fR,chi,[0,0],'linear');

idxR   = (fR>0.02) & (fR<0.1);
W_azi = mean(W(:,idxR),2,'omitnan')';

W_azi_cs = (W_azi-mean(W_azi))./std(W_azi);

A1 = lsqcurvefit(fitfun,[-1,3,0.9,4,0.5,0.05,1.5*pi],chi/180*pi,W_azi_cs);
mu_degrees   = A1(3)*180/pi;
fwhm_degrees = 2*sqrt(2*log(2)/A1(4))/pi*180;

W_s   = unwrapPolar(fX,fY,M_s,fR,chi,[0,0],'linear');
W_sbkg        = mean(W_s,1,'omitnan')';

wedge_angle = 15;
chi_idx_fiber_norm = abs(chi-mu_degrees)<=wedge_angle;
W_sfiber_norm      = mean(W_s(chi_idx_fiber_norm,:),1,'omitnan')';

chi_idx_fiber = abs(chi-90-mu_degrees)<=wedge_angle;
W_sfiber = mean(W_s(chi_idx_fiber,:),1,'omitnan')';

%% Extract line profile oriented normal to predominant fiber orientation

lp_x0 = 93;
lp_y0 = 52.2;
dir1 = RotM(mu_degrees)*[1;0];
dir2 = RotM(mu_degrees)*[0;1];

lp_xc = lp_x0 + 30*dir1(1);
lp_yc = lp_y0 + 30*dir1(2);

xqv = linspace(-30,30,300);
yqv = linspace(-10,+10,100);

[Xq,Yq]=meshgrid(xqv,yqv);
XYqr = RotM(mu_degrees)*[Xq(:).';Yq(:).'];


% use interpolation to extract area for line profile and integrate for
% better S/N
ROI_I_sel = reshape(interp2(xv_sel,yv_sel,I_sel,XYqr(1,:)+lp_xc,XYqr(2,:)+lp_yc),100,300);

% generate vertices and faces for a patch to indicate selection
Verts_BoxBLM = (RotM(mu_degrees)*[[xqv(1);yqv(1)],[xqv(end);yqv(1)],[xqv(end);yqv(end)],[xqv(1);yqv(end)]]+[lp_xc;lp_yc]).'; % [nm]
Faces_BoxBLM = [1,2,3,4];

%%  Extract line profiles from smoothed FFT magnitude 

lp_fft_x0 = 0;
lp_fft_y0 = 0;
lp_fft_x1 = max(freq);
lp_fft_y1 = max(freq);

lp_fft_1 = RotM(mu_degrees-45)*[[lp_fft_x0,lp_fft_x1];[lp_fft_y0,lp_fft_y1]];
lp_fft_2 = RotM(90)*lp_fft_1;

[cx11,cy12,c1] = improfile(freq,freq,M_s,lp_fft_1(1,:),lp_fft_1(2,:),1000,'bilinear');
[cx21,cy22,c2] = improfile(freq,freq,M_s,lp_fft_2(1,:),lp_fft_2(2,:),1000,'bilinear');

rc1 = sqrt(cx11.^2+cy12.^2);
rc2 = sqrt(cx21.^2+cy22.^2);

c1 = c1(~isnan(c1));
c2 = c2(~isnan(c2));
rc1 = rc1(~isnan(c1));
rc2 = rc2(~isnan(c2));

c2 = c2(1:length(c1));
rc2 = rc2(1:length(c1));
%% Radial avg
close all

dr = max(rc1)/(length(rc1)-1);
edges   = -dr/2:dr:max(rc1)+dr/2;
r_int = edges(1:end-1)+dr/2;

[~,M_rint] = radial_avg(fX,fY,M_s,edges);

%% Remove estimated background from smoothed FFT magnitude
% create square grid over frequency domain 
[fY,fX] = ndgrid(freq,freq);

% calculate distance of gridpoints from origin to use as query points for
% interpolation; note query points are a N-by-2 vector
fRq = sqrt((fX(:)).^2+(fY(:)).^2);

% calculate background magnitude from off-diagonal profile by creating
% interpolated surface of revolution; note that output is a N-by-1 vector!
M_Bkg = interp1(rc2,c2,fRq);
M_Bkg2 = interp1(r_int,M_rint,fRq);
M_Bkg3 = interp1(fR,W_sbkg,fRq);

% reshape vector to array
M_Bkg = reshape(M_Bkg,N,N);
M_Bkg2 = reshape(M_Bkg2,N,N);
M_Bkg3 = reshape(M_Bkg3,N,N);

% subtract background from smoothed magnitude 
M_corr  = M_s-M_Bkg;
M_corr2 = M_s-M_Bkg2;
M_corr3 = M_s-M_Bkg3;

% shift corrected magnitude to make sure values are >=0 
M_c = M_corr-min(M_corr,[],'all');
M_c2 = M_corr2-min(M_corr2,[],'all');
M_c3 = M_corr3-min(M_corr3,[],'all');

%% apply circular Butterworth filter with radius R
[ri, ci] = find(M_corr == max(M_corr,[],'all'));

R   = 0.06;
dirv = (RotM(mu_degrees)*[1;0]).'; % direction vector of cyan line in FFT

if R > 0.05
    Ori = freq([ci,ri])-(R-0.05)*dirv; % to avoid overlapping circles shift origin along diagonal
else
    Ori = freq([ci,ri]);
end
fD1 = sqrt((fX-Ori(1)).^2+(fY-Ori(2)).^2);
fD2 = sqrt((fX+Ori(1)).^2+(fY+Ori(2)).^2);

order = 4;

H = butter_fun(fD1,R,order)+butter_fun(fD2,R,order);

F_f = H.*F;
I_f = ifft2(ifftshift(F_f));

%% apply Butterworth lowpass filter

fD1 = sqrt((fX).^2+(fY).^2);

ubf = norm(freq([ci,ri]))+norm((R-0.05)*dirv)+R; 
H_lp = butter_fun(fD1,ubf,order);
F_lp = H_lp.*F;
I_lp = ifft2(ifftshift(F_lp));

%% apply BW bandpass filter
H_bp = butter_fun(fD1,0.02,order).*(1-butter_fun(fD1,0.3,order));
F_bp = H_bp.*F;
I_bp = ifft2(ifftshift(F_bp));

%% extract line profiles in filtered images 
ROI_I_f  = reshape(interp2(xv_sel,yv_sel,abs(I_f),XYqr(1,:)+lp_xc,XYqr(2,:)+lp_yc),100,300);
ROI_I_lp = reshape(interp2(xv_sel,yv_sel,abs(I_lp),XYqr(1,:)+lp_xc,XYqr(2,:)+lp_yc),100,300);
ROI_I_bp = reshape(interp2(xv_sel,yv_sel,I_bp,XYqr(1,:)+lp_xc,XYqr(2,:)+lp_yc),100,300);

%% visualize
close all
ctr = 1;
panel = 'ABCDEFGHIKLMNOPQRST';

f1 = figure('Position',[150 250 1850 1450]);
tile1 = tiledlayout(3,4);

% Panel A
ctr_img=ctr;
ax{ctr} = nexttile;
imagesc(xv,yv,I);hold on;
axis equal tight
patch('Faces',Faces_BoxA,'Vertices',Verts_BoxA,'FaceColor','none','EdgeColor','w','LineWidth',1);
ax{ctr}.TickDir = 'out';
xlabel('x [nm]');
ylabel('y [nm]');
title([panel(ctr),'. original image of stylus composite']);

% Panel B
ctr = ctr+1;
ctr_roi = ctr;
ax{ctr} = nexttile;
imagesc(xv_sel,yv_sel,I_sel);hold on;
patch('Faces',Faces_BoxBLM,'Vertices',Verts_BoxBLM,'FaceColor','none','EdgeColor','w','LineWidth',1);

ax{ctr}.TickDir = 'out';
xlabel('x [nm]');
ylabel('y [nm]');
axis equal tight
colormap(gray)
title([panel(ctr),'. closeup of area indicated in ',panel(ctr_img)]);

% Panel C
ctr = ctr+1;
ctr_F = ctr;
ax{ctr} = nexttile;
imagesc(freq,freq,log(M));hold on;
ax{ctr}.TickDir = 'out';
axis equal tight
xlabel('frequency [nm^{-1}]');
ylabel('frequency [nm^{-1}]');
colormap(gray)
title([panel(ctr),'. FFT of (',panel(ctr_roi),').']);

% Panel D
ctr = ctr+1;
ax{ctr} = nexttile;
plot(chi,W_azi_cs,'.');hold on;

plot(chi,fitfun(A1,chi*pi/180),'r','LineWidth',1);
plot(mu_degrees*[1,1],fitfun(A1,A1(3))*[1,1.5],'r','LineWidth',1);
plot(mu_degrees+fwhm_degrees*[-0.5,0.5],fitfun([A1(1:4),0,1,1],(mu_degrees+fwhm_degrees*[-0.5,0.5])*pi/180),'r','LineWidth',1);
text(mu_degrees,fitfun(A1,A1(3))*1.55,num2str(mu_degrees,'%2.2f'),'HorizontalAlignment','center');
text(mu_degrees,fitfun(A1,(mu_degrees-fwhm_degrees/2)*pi/180)*0.9,num2str(fwhm_degrees,'%2.1f'),'HorizontalAlignment','center');

ax{ctr}.TickDir = 'out';
ax{ctr}.XLim = [0,180];
xlabel('\theta [^\circ]');
ylabel('accumulated intensity');
title([panel(ctr),'. Orientation Distribution']);

% Panel E
ctr = ctr+1;
ctr_Fs = ctr; 
ax{ctr} = nexttile;
imagesc(freq,freq,log(M_s));hold on;
axis equal tight
rotate(plot(lp_fft_1(1,:),lp_fft_1(2,:),'c'),[0,0,1],wedge_angle,[0,0,0]);
rotate(plot(lp_fft_1(1,:),lp_fft_1(2,:),'c'),[0,0,1],-wedge_angle,[0,0,0]);
rotate(plot(lp_fft_2(1,:),lp_fft_2(2,:),'r'),[0,0,1],wedge_angle,[0,0,0]);
rotate(plot(lp_fft_2(1,:),lp_fft_2(2,:),'r'),[0,0,1],-wedge_angle,[0,0,0]);

colormap(gray)
ax{ctr}.TickDir = 'out';
ax{ctr}.XLim = [min(freq),max(freq)];
ax{ctr}.YLim = [min(freq),max(freq)];
xlabel('frequency [nm^{-1}]');
ylabel('frequency [nm^{-1}]');
title([panel(ctr),'. Smoothed FFT of (',panel(ctr_roi),').']);

% Panel F
ctr = ctr+1;
ctr_lp = ctr;
ax{ctr} = nexttile;

loglog(fR,W_sfiber_norm,'c.'); hold on;
loglog(fR,W_sfiber,'r.');
loglog(fR,W_sbkg,'k.'); 

xlabel('frequency [nm^{-1}]');
ylabel('intensity [a.u.]');
ax{ctr}.TickDir = 'out';
ax{ctr}.XLim = [1e-2,5];
title([panel(ctr),'. Wedge integrals indicated in (',panel(ctr_Fs),').']);
legend('fibre normal','fibre axis','radial avg')

% Panel G
ctr = ctr+1;
ax{ctr} = nexttile;

plot(fR,W_sfiber_norm-W_sfiber,'k.'); hold on;
plot(fR,W_sfiber_norm-W_sbkg,'r.');
ax{ctr}.TickDir = 'out';
ax{ctr}.XTick = [1/16,1/10,1/7,1/6,1/5,1/4,1/3,1/2];
ax{ctr}.XTickLabels = {'16^{-1}','10^{-1}','7^{-1}','6^{-1}','5^{-1}','4^{-1}','3^{-1}','2^{-1}'}
ax{ctr}.XLim = [0.01,0.5];
ax{ctr}.YLim = [0,3e5];
xlabel('frequency [nm^{-1}]');
ylabel('intensity [a.u.]');
title([panel(ctr),'. Difference of wedge integrals in (',panel(ctr_lp),').']);
legend('fibre normal - fibre axis','fibre normal - radial average')
grid on

% Panel H
ctr = ctr+1;
ctr_bkg = ctr;
ax{ctr} = nexttile;
imagesc(freq,freq,log(M_Bkg3));hold on;
ax{ctr}.TickDir = 'out';
axis equal tight
xlabel('frequency [nm^{-1}]');
ylabel('frequency [nm^{-1}]');
colormap(gray)
title({[panel(ctr),'. Background Estimate:'],['Surface of Revolution of black profile in (',panel(ctr_lp),').']});

% Panel I
ctr = ctr+1;
ctr_cor = ctr;
ax{ctr} = nexttile;
imagesc(freq,freq,M_corr3);hold on
ax{ctr}.TickDir = 'out';
axis equal tight
rotate(plot(lp_fft_1(1,:),lp_fft_1(2,:),'c'),[0,0,1],wedge_angle,[0,0,0]);
rotate(plot(lp_fft_1(1,:),lp_fft_1(2,:),'c'),[0,0,1],-wedge_angle,[0,0,0]);
rotate(plot(lp_fft_2(1,:),lp_fft_2(2,:),'r'),[0,0,1],wedge_angle,[0,0,0]);
rotate(plot(lp_fft_2(1,:),lp_fft_2(2,:),'r'),[0,0,1],-wedge_angle,[0,0,0]);
xlim([min(freq),max(freq)]);
ylim([min(freq),max(freq)]);
xlabel('frequency [nm^{-1}]');
ylabel('frequency [nm^{-1}]');
title([panel(ctr),'. Background-corrected FFT: (',panel(ctr_Fs),') - (',panel(ctr_bkg),').']);

%  Panel K 
ctr = ctr+1;
ax{ctr} = nexttile;
imagesc(freq,freq,M_corr3); hold on;
xlabel('frequency [nm^{-1}]');
ylabel('frequency [nm^{-1}]');
axis equal tight
rotate(plot(lp_fft_1(1,:),lp_fft_1(2,:),'c'),[0,0,1],wedge_angle,[0,0,0]);
rotate(plot(lp_fft_1(1,:),lp_fft_1(2,:),'c'),[0,0,1],-wedge_angle,[0,0,0]);
rotate(plot(lp_fft_2(1,:),lp_fft_2(2,:),'r'),[0,0,1],wedge_angle,[0,0,0]);
rotate(plot(lp_fft_2(1,:),lp_fft_2(2,:),'r'),[0,0,1],-wedge_angle,[0,0,0]);
% indicate ideal filter for spot filter
viscircles([Ori;-Ori],[R,R],'Color','y','LineStyle',':','LineWidth',1,'EnhanceVisibility',false);

% indicate ideal filter for lowpass
viscircles([0,0],ubf,'Color','m','LineStyle',':','LineWidth',1,'EnhanceVisibility',false);
ax{ctr}.TickDir = 'out';
ax{ctr}.XLim = [-0.4,0.4];
ax{ctr}.YLim = [-0.4,0.4];
xlabel('frequency [nm^{-1}]');
ylabel('frequency [nm^{-1}]');
title([panel(ctr),'. Close-up of (',panel(ctr_cor),').']);

% Panel L
ctr = ctr+1;
ctr_bw = ctr;
ax{ctr} = nexttile;
imagesc(xv_sel,yv_sel,abs(I_f));hold on;
patch('Faces',Faces_BoxBLM,'Vertices',Verts_BoxBLM,'FaceColor','none','EdgeColor','w','LineWidth',1);
xlabel('x [nm]');
ylabel('y [nm]');
ax{ctr}.TickDir = 'out';
axis equal tight
title({[panel(ctr),'. Inverse FFT after appl. of Butterworth'],['filter to (',panel(ctr_F),'); yellow circles in (K).']});

% Panel M
ctr = ctr+1;
ctr_bp = ctr;
ax{ctr} = nexttile;
imagesc(xv_sel,yv_sel,abs(I_lp));hold on;

patch('Faces',Faces_BoxBLM,'Vertices',Verts_BoxBLM,'FaceColor','none','EdgeColor','w','LineWidth',1);
xlabel('x [nm]');
ylabel('y [nm]');
ax{ctr}.TickDir = 'out';
axis equal tight
title({[panel(ctr),'. Inverse FFT after appl. of Butterworth'],[' lowpass filter to (',panel(ctr_F),'); magenta circle in (K).']});

%% Figure 2

original = mean(ROI_I_sel,1);
original = original - mean(original);
masked   = mean(ROI_I_f,1);
masked   = masked - mean(masked);
lowpass  = mean(ROI_I_lp,1);
lowpass  = lowpass - mean(lowpass);

panel = 'ABCDEFGHIKLMNOPQRSTU';
f2 = figure('Position',[1124 841 392 946]);
tile2 = tiledlayout(5,1);

ctr = ctr+1;
ax{ctr} = nexttile;
imagesc(xqv,yqv,ROI_I_sel);
daspect([1 1 1]);
ax{ctr}.XTick =[];
ax{ctr}.TickDir = 'out';
ylabel('y [nm]');
title([panel(ctr),'. ROI in (',panel(ctr_roi),').']);

ctr = ctr+1;
ax{ctr} = nexttile;
imagesc(xqv,yqv,ROI_I_f);
daspect([1 1 1]);
ax{ctr}.XTick =[];
ax{ctr}.TickDir = 'out';
ylabel('y [nm]');
title([panel(ctr),'. ROI in (',panel(ctr_bw),').']);

ctr = ctr+1;
ax{ctr} = nexttile;
imagesc(xqv,yqv,ROI_I_lp);
daspect([1 1 1]);
ax{ctr}.TickDir = 'out';
ax{ctr}.XMinorTick = 'on';
colormap gray
xlabel('x [nm]');
ylabel('y [nm]');
title([panel(ctr),'. ROI in (',panel(ctr_bp),').']);

ctr = ctr+1;
ax{ctr} = nexttile([2,1]);

plot(xqv,original,'LineWidth',1); hold on;
plot(xqv,masked,'LineWidth',1);
plot(xqv,lowpass,'LineWidth',1);
x{ctr}.TickDir = 'out';
ax{ctr}.XMinorTick = 'on';
ax{ctr}.YTick =[];
colormap gray
xlabel('x [nm]');
ylabel('I [a.u.]');
legend('original','circ. bw mask','lowpass bw filter')
title([panel(ctr),'. Integrated Profiles in (',panel(ctr_roi),'), (',panel(ctr_bw),'), and (',panel(ctr_bp),').']);
%%
if save_figure
    saveas(f1,[pn_out,fn_out{1}],'epsc');
    disp(['Saved to file:',pn_out,fn_out{1}] );
    saveas(f2,[pn_out,fn_out{2}],'epsc');
    disp(['Saved to file:',pn_out,fn_out{2}] );
else
    disp('Figures not saved!');
end

%%

function unwrapped = unwrapPolar(X, Y, I, rho, chi, offset,method)
% rho is a vector of radial positions
% chi is a vector of azimuthal angles
% offset is the displacement of the beam center from the actual beam center
       
    % define ND grid of query rays
    [Rhoq,Chiq] = meshgrid(rho,chi);
    [Xq,Yq] = pol2cart(Chiq/180*pi,Rhoq);     
    % interpolate
    unwrapped = interp2(X,Y,I,Xq+offset(1),Yq+offset(2),method,NaN);
end

function warped = warpCart(Rho, Chi, I, x, y, offset)
% x is a vector of x positions
% y is a vector of y positions
% offset is the displacement of the beam center from the actual beam center
       
    % define ND grid of query rays
    [Xq,Yq] = meshgrid(x,y);
    [Chiq,Rhoq] = cart2pol(Xq,Yq);     
    % interpolate
    warped = interp2(Rho,Chi,I,Rhoq,(Chiq+pi)/pi*180);
end

function rS = radial_sum(X,Y,A,m)
    [Theta, Rho] = cart2pol(X,Y);
    [~,edges,bins] = histcounts(Rho,m);
    centers = cumsum(mean(edges([1,2]))+diff(edges));
    S = accumarray(bins(:),A(:));
end

function [edges,rAvg] = radial_avg(X,Y,A,m)
    [Theta, Rho] = cart2pol(X,Y);
    [~,edges,bins] = histcounts(Rho,m);
    centers = cumsum(mean(edges([1,2]))+diff(edges));
    rAvg = accumarray(reshape(bins(bins>0),[],1),reshape(A(bins>0),[],1),[],@(x) mean(x,'omitnan'));
end
