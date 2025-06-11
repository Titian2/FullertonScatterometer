%% --------------------------------------------------------------------------
%  Scattering-map simulator – Gaussian vs. flat-top illumination intercept
% --------------------------------------------------------------------------
clear; close all; clc;

%% ---------------- USER PARAMETERS (unchanged) -----------------------------
sample_diam_mm  = 25.4;
beam_FWHM_mm    = [5 5];
beam_center_mm  = [0 0];
tilt_deg        = 10;
num_pixels      = 1024;

num_scatterers  = 1000;
ROI_radii_mm    = linspace(2.5,3.5,6);
max_iter        = 1000;

inital_scatter_diam = 0.005;                % mm
growth_mm           = inital_scatter_diam*0.001;
bright_step         = 0.10;
dark_level          = 0;

beam_threshold = 0.10;                      % ≥10 % of peak

%% ---------------- GRID & STATIC MAPS (unchanged) --------------------------
px_per_mm  = num_pixels / sample_diam_mm;
[Xpx,Ypx]  = meshgrid(1:num_pixels,1:num_pixels);
x_mm       = (Xpx-(num_pixels+1)/2)/px_per_mm;
y_mm       = (Ypx-(num_pixels+1)/2)/px_per_mm;
sample_mask = (x_mm.^2 + y_mm.^2) <= (sample_diam_mm/2)^2;

d0        = 100;
tilt_map  = (d0./(d0 + x_mm*sind(tilt_deg))).^2 .* sample_mask;

sig_x = beam_FWHM_mm(1)/(2*sqrt(2*log(2)));
sig_y = beam_FWHM_mm(2)/(2*sqrt(2*log(2)));
beam_map = exp( -((x_mm-beam_center_mm(1)).^2)/(2*sig_x^2) ...
                -((y_mm-beam_center_mm(2)).^2)/(2*sig_y^2) );
beam_map(~sample_mask) = 0;

illum_gauss = beam_map .* tilt_map;                 % true profile
beam_mask   = illum_gauss >= beam_threshold*max(illum_gauss(:));

illum_flat  = double(beam_mask);      % strictly uniform intensity
Ptot_flat   = sum(illum_flat(:));     % its total power

Ptot_gauss  = sum(illum_gauss(beam_mask));          % normalisations

%% -------- INITIAL SCATTERERS (generated only inside beam) -----------------
rng default
beam_lin  = find(beam_mask);
pick_inds = beam_lin(randi(numel(beam_lin),1,num_scatterers));
[py, px]  = ind2sub([num_pixels num_pixels],pick_inds);
sc_cx = (px-(num_pixels+1)/2)/px_per_mm;
sc_cy = (py-(num_pixels+1)/2)/px_per_mm;
sc_r  = inital_scatter_diam*ones(1,num_scatterers);
sc_amp= rand(1,num_scatterers)*20;

scatterMap = zeros(num_pixels,'single');

%% --------------------------------------------------------------------------
%  Scattering-map simulator – Gaussian vs. flat-top illumination intercept
% --------------------------------------------------------------------------
clear; close all; clc;

%% ---------------- USER PARAMETERS (unchanged) -----------------------------
sample_diam_mm  = 25.4;
beam_FWHM_mm    = [5 5];
beam_center_mm  = [0 0];
tilt_deg        = 10;
num_pixels      = 1024;

num_scatterers  = 1000;
ROI_radii_mm    = linspace(2.5,3.5,6);
max_iter        = 1000;

inital_scatter_diam = 0.005;                % mm
growth_mm           = inital_scatter_diam*0.001;
bright_step         = 0.10;
dark_level          = 0;

beam_threshold = 0.10;                      % ≥10 % of peak

%% ---------------- GRID & STATIC MAPS (unchanged) --------------------------
px_per_mm  = num_pixels / sample_diam_mm;
[Xpx,Ypx]  = meshgrid(1:num_pixels,1:num_pixels);
x_mm       = (Xpx-(num_pixels+1)/2)/px_per_mm;
y_mm       = (Ypx-(num_pixels+1)/2)/px_per_mm;
sample_mask = (x_mm.^2 + y_mm.^2) <= (sample_diam_mm/2)^2;

d0        = 100;
tilt_map  = (d0./(d0 + x_mm*sind(tilt_deg))).^2 .* sample_mask;

sig_x = beam_FWHM_mm(1)/(2*sqrt(2*log(2)));
sig_y = beam_FWHM_mm(2)/(2*sqrt(2*log(2)));
beam_map = exp( -((x_mm-beam_center_mm(1)).^2)/(2*sig_x^2) ...
                -((y_mm-beam_center_mm(2)).^2)/(2*sig_y^2) );
beam_map(~sample_mask) = 0;

illum_gauss = beam_map .* tilt_map;                 % true profile


beam_mask   = illum_gauss >= beam_threshold*max(illum_gauss(:));

illum_flat  = double(beam_mask);      % strictly uniform intensity


%% -------- INITIAL SCATTERERS (generated only inside beam) -----------------
rng default
beam_lin  = find(beam_mask);
pick_inds = beam_lin(randi(numel(beam_lin),1,num_scatterers));
[py, px]  = ind2sub([num_pixels num_pixels],pick_inds);
sc_cx = (px-(num_pixels+1)/2)/px_per_mm;
sc_cy = (py-(num_pixels+1)/2)/px_per_mm;
sc_r  = inital_scatter_diam*ones(1,num_scatterers);
sc_amp= rand(1,num_scatterers)*20;

scatterMap = zeros(num_pixels,'single');


for k = 1:num_scatterers
    m = ((x_mm-sc_cx(k)).^2 + (y_mm-sc_cy(k)).^2)<= (sc_r(k)*px_per_mm)^2;
    scatterMap(m) = scatterMap(m) + sc_amp(k);
end

% ---- create two IDENTICAL copies of the initial scatter state -------------
scatterG = scatterMap;           % will be multiplied by illum_gauss
scatterF = scatterMap;           % independent copy for flat-top


% ---- power totals ----------------------------------------------------------
Ptot_gauss = sum(illum_gauss(:));            % Gaussian normalisation
Ptot_flat  = sum(illum_flat(:));             % just #pixels inside beam


%% ---------------- ROI MASKS (unchanged) -----------------------------------
ROI_masks = false(num_pixels,num_pixels,numel(ROI_radii_mm));
for i = 1:numel(ROI_radii_mm)
    ROI_masks(:,:,i) = ((x_mm-beam_center_mm(1)).^2 + ...
                        (y_mm-beam_center_mm(2)).^2) <= ROI_radii_mm(i)^2;
end


%% ---------------- MAIN LOOP -----------------------------------------------
yG = nan(max_iter,1,'single');    % Gaussian intercepts
yF = nan(max_iter,1,'single');    % Flat-top intercepts

for it = 1:max_iter
    %--- grow & brighten BOTH scatter maps identically ----------------------
    for k = 1:num_scatterers
        sc_r(k) = sc_r(k) + growth_mm;
        mask_k  = ((x_mm-sc_cx(k)).^2 + (y_mm-sc_cy(k)).^2) ...
                  <= (sc_r(k)*px_per_mm)^2;
        scatterG(mask_k) = scatterG(mask_k) + bright_step;
        scatterF(mask_k) = scatterF(mask_k) + bright_step;
    end
    
    %--- physics frames -----------------------------------------------------
    imgG = illum_gauss .* scatterG + dark_level*(sample_mask & ~beam_mask);
    imgF = illum_flat  .* scatterF + dark_level*(sample_mask & ~beam_mask);
    
    %--- live image shows Gaussian-weighted scatter (or use scatterG) -------
    set(hImg,'CData',scatterG);     % purely cosmetic
    
    %--- intercepts (Gaussian vs. flat-top) --------------------------------
    for i = 1:numel(ROI_radii_mm)
        m        = ROI_masks(:,:,i);
        sumsG(i) = sum(imgG(m));    sumsF(i) = sum(imgF(m));
        areas(i) = sum(m(:));
    end
    
    coeffG   = polyfit(double(areas),double(sumsG),1);
    yG(it)   = coeffG(2) / Ptot_gauss;
    coeffF   = polyfit(double(areas),double(sumsF),1);
    yF(it)   = coeffF(2) / Ptot_flat;
    
    set(hPlotG,'XData',1:it,'YData',yG(1:it));
    set(hPlotF,'XData',1:it,'YData',yF(1:it));
    drawnow;
end




for k = 1:num_scatterers
    m = ((x_mm-sc_cx(k)).^2 + (y_mm-sc_cy(k)).^2)<= (sc_r(k)*px_per_mm)^2;
    scatterMap(m) = scatterMap(m) + sc_amp(k);
end

%% ---------------- ROI MASKS (unchanged) -----------------------------------
ROI_masks = false(num_pixels,num_pixels,numel(ROI_radii_mm));
for i = 1:numel(ROI_radii_mm)
    ROI_masks(:,:,i) = ((x_mm-beam_center_mm(1)).^2 + ...
                        (y_mm-beam_center_mm(2)).^2) <= ROI_radii_mm(i)^2;
end

%% ---------------- FIGURE SET-UP (only intercept plot changed) -------------
f = figure('Color','w','Units','normalized','Position',[0.05 0.1 0.9 0.75]);

hAxImg = axes(f,'Position',[0.08 0.25 0.35 0.65]);
hImg   = imagesc(scatterMap,'Parent',hAxImg); axis(hAxImg,'image','off');
colormap(hAxImg,'bone'); caxis(hAxImg,[0 max(scatterMap(:))+eps]); hold on;

theta = linspace(0,2*pi,360);  % ROI & rim (unchanged)
for i = 1:numel(ROI_radii_mm)
    r_px = ROI_radii_mm(i)*px_per_mm;
    plot(hAxImg,(num_pixels+1)/2+r_px*cos(theta), ...
                 (num_pixels+1)/2+r_px*sin(theta),'m--','LineWidth',1.2);
end
r_px_disk = (sample_diam_mm/2)*px_per_mm;
plot(hAxImg,(num_pixels+1)/2+r_px_disk*cos(theta), ...
             (num_pixels+1)/2+r_px_disk*sin(theta),'w--','LineWidth',1.1);
title(hAxImg,'Scattering Simulation');

pos = get(hAxImg,'Position');                                % beam profiles

% Reduce the height of the vertical axes (hAxY)
y_height = pos(4) * 0.7;   % 70% of original height
y_bottom = pos(2) + (pos(4) - y_height)/2;
hAxY = axes(f,'Position',[pos(1)-0.09 y_bottom 0.08 y_height]);
plot(hAxY,illum_gauss(:,num_pixels/2),1:num_pixels,'k');
set(hAxY,'YDir','reverse','Box','off');

% Decrease the length of the horizontal axes (hAxX)
x_width = pos(3) * 0.7;    % 70% of original width
x_left = pos(1) + (pos(3) - x_width)/2;
hAxX = axes(f,'Position',[x_left pos(2)-5.50 x_width 0.08]);
plot(hAxX,1:num_pixels,illum_gauss(num_pixels/2,:),'k'); set(hAxX,'Box','off');

set(hAxY, ...
    'Position', [0.017  0.24  0.05  0.66], ...
    'YLim',    [0 1000], ...        % show rows   0–1000 only
    'xlabel',  "Y- Projection");


set(hAxX, ...
    'Position', [0.08  0.15  0.35  0.08], ...
    'XLim',    [0 1000], ...        % show columns 0–1000 only
    'xlabel' , "X- Projection");

linkaxes([hAxImg hAxX],'x'); linkaxes([hAxImg hAxY],'y');

% -------- new dual-curve intercept plot ------------------------------------
hAxInt = axes(f,'Position',[0.55 0.25 0.40 0.65]); hold on;
hPlotG = plot(hAxInt,nan,nan,'b-','LineWidth',1.4,'DisplayName','Gaussian');
hPlotF = plot(hAxInt,nan,nan,'r--','LineWidth',1.4,'DisplayName','Flat-top');
grid(hAxInt,'on'); xlabel('Iteration'); ylabel('Intercept / P_{tot}');
title(hAxInt,'Blue = Gaussian    Red = Flat-top');
legend(hAxInt,'Location','northwest');

%% ---------------- MAIN LOOP -----------------------------------------------
yG = nan(max_iter,1,'single');    % Gaussian intercepts
yF = nan(max_iter,1,'single');    % Flat-top intercepts

for it = 1:max_iter
    % --- grow & brighten defects
    for k = 1:num_scatterers
        sc_r(k) = sc_r(k) + growth_mm;
        m = ((x_mm-sc_cx(k)).^2 + (y_mm-sc_cy(k)).^2)<= (sc_r(k)*px_per_mm)^2;
        scatterMap(m) = scatterMap(m) + bright_step;
    end
    
    % --- frames for Gaussian vs. flat-top weighting
    imgG = illum_gauss .* scatterMap + dark_level*(sample_mask.*~beam_mask);
    imgF = illum_flat  .* scatterMap + dark_level*(sample_mask & ~beam_mask);
    
    set(hImg,'CData',scatterMap);                 % live view
    
    % --- intercepts for each ROI set
    sumsG = zeros(numel(ROI_radii_mm),1,'single');
    sumsF = zeros(size(sumsG),'single');
    areas = zeros(size(sumsG),'single');
    for i = 1:numel(ROI_radii_mm)
        m = ROI_masks(:,:,i);
        sumsG(i) = sum(imgG(m));  sumsF(i) = sum(imgF(m));
        areas(i) = sum(m(:));
    end
    
    % linear fits → y-intercepts → normalised
    coeffG      = polyfit(double(areas), double(sumsG), 1);
    yG(it)      = coeffG(2) / Ptot_gauss;      % intercept / power
    
    coeffF      = polyfit(double(areas), double(sumsF), 1);
    yF(it)      = coeffF(2) / Ptot_flat;       % intercept / power
    
    % update curves
    set(hPlotG,'XData',1:it,'YData',yG(1:it));
    set(hPlotF,'XData',1:it,'YData',yF(1:it));
    drawnow;
end