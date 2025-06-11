%% --------------------------------------------------------------------------
%  Scattering-map simulator – Gaussian vs. flat-top illumination intercept
% --------------------------------------------------------------------------
clear; close all; clc;

%% 1. USER PARAMETERS --------------------------------------------------------
sample_diam_mm   = 25.4;                   % disk diameter (mm)
beam_FWHM_mm     = [5 5];                  % Gaussian FWHM (x,y) in mm
beam_center_mm   = [0 0];                  % beam offset wrt disk centre
tilt_deg         = 10;                     % left→right 1/d^2 fall-off (deg)
num_pixels       = 1024;                   % simulation grid N×N

num_scatterers   = 1000;                   % number of scatterers
initial_diam_mm  = 0.005;                  % starting scatter diameter (mm)
growth_mm        = initial_diam_mm * 0.001;% radius growth per iteration (mm)
bright_step      = 0.10;                   % brightness boost per iteration
max_iter         = 1000;                   % total iterations

beam_threshold   = 0.10;                   % ≥10% of peak → inside beam
dark_level       = 0;                      % background outside beam

ROI_radii_mm     = linspace(2.5,3.5,6);    % ROI radii in mm

%% 2. GRID & STATIC MAPS -----------------------------------------------------
px_per_mm  = num_pixels / sample_diam_mm;
[Xpx,Ypx]  = meshgrid(1:num_pixels, 1:num_pixels);
x_mm       = (Xpx - (num_pixels+1)/2) / px_per_mm;
y_mm       = (Ypx - (num_pixels+1)/2) / px_per_mm;

disk_mask = (x_mm.^2 + y_mm.^2) <= (sample_diam_mm/2)^2;

% Tilt map
d0       = 1000;  % nominal camera distance (mm)
tilt_map = (d0 ./ (d0 + x_mm * sind(tilt_deg))).^2 .* disk_mask;

% Gaussian beam map
sig_x    = beam_FWHM_mm(1) / (2*sqrt(2*log(2)));
sig_y    = beam_FWHM_mm(2) / (2*sqrt(2*log(2)));
beam_map = exp( -((x_mm - beam_center_mm(1)).^2)/(2*sig_x^2) ...
                -((y_mm - beam_center_mm(2)).^2)/(2*sig_y^2) );

beam_map(~disk_mask) = 0.1;

% Illumination profiles
illum_gauss = beam_map .* tilt_map;                  % true Gaussian × tilt

beam_mask   = illum_gauss >= beam_threshold * max(illum_gauss(:));
illum_flat  = double(beam_mask);                     % uniform inside beam

Ptot_gauss  = sum(illum_gauss(:));                   % normalisation
Ptot_flat   = sum(illum_flat(:));

%% 3. INITIAL DEFECT MAP (inside beam only) ---------------------------------
rng default
beam_lin   = find(beam_mask);
pick_inds  = beam_lin(randi(numel(beam_lin),1,num_scatterers));
[py, px]   = ind2sub([num_pixels num_pixels], pick_inds);
sc_cx      = (px - (num_pixels+1)/2) / px_per_mm;
sc_cy      = (py - (num_pixels+1)/2) / px_per_mm;
sc_r       = initial_diam_mm * ones(1,num_scatterers);
sc_amp     = rand(1,num_scatterers) * 20;

scatterG   = zeros(num_pixels,'single');  % for Gaussian branch
scatterF   = zeros(num_pixels,'single');  % for flat-top branch
for k = 1:num_scatterers
    mask = ((x_mm - sc_cx(k)).^2 + (y_mm - sc_cy(k)).^2) ...
           <= (sc_r(k) * px_per_mm)^2;
    scatterG(mask) = scatterG(mask) + sc_amp(k);
    scatterF(mask) = scatterF(mask) + sc_amp(k);
end

%% 4. ROI MASKS --------------------------------------------------------------
ROI_masks = false(num_pixels, num_pixels, numel(ROI_radii_mm));
for i = 1:numel(ROI_radii_mm)
    ROI_masks(:,:,i) = ((x_mm - beam_center_mm(1)).^2 + ...
                        (y_mm - beam_center_mm(2)).^2) ...
                        <= (ROI_radii_mm(i))^2;
end

%% 5. FIGURE LAYOUT ----------------------------------------------------------
f = figure('Color','w','Units','normalized','Position',[0.05 0.1 0.9 0.75]);

% Main scatter image (Gaussian branch displayed)
hAxImg = axes(f,'Position',[0.08 0.25 0.35 0.65]);
hImg   = imagesc(scatterG,'Parent',hAxImg);
axis(hAxImg,'image','off');
colormap(hAxImg,'bone');
caxis(hAxImg,[0 max(scatterG(:))+eps]);
hold(hAxImg,'on');

% ROI rings and disk rim
theta     = linspace(0,2*pi,360);
for i = 1:numel(ROI_radii_mm)
    r_px = ROI_radii_mm(i) * px_per_mm;
    plot(hAxImg, (num_pixels+1)/2 + r_px*cos(theta), ...
                 (num_pixels+1)/2 + r_px*sin(theta), 'm--', 'LineWidth',1.2);
end
r_px_disk = (sample_diam_mm/2) * px_per_mm;
plot(hAxImg, (num_pixels+1)/2 + r_px_disk*cos(theta), ...
             (num_pixels+1)/2 + r_px_disk*sin(theta), 'w--', 'LineWidth',1.1);
title(hAxImg,'Simulated Scatter');

% Y-profile axis (manual position, cropped 0–1000)
hAxY = axes(f,'Position',[0.017 0.24 0.05 0.66]);
plot(hAxY, illum_gauss(:,num_pixels/2), 1:num_pixels, 'k');
set(hAxY, 'YDir','reverse', 'Box','off', 'YLim', [0 1000]);

% X-profile axis (manual position, cropped 0–1000)
hAxX = axes(f,'Position',[0.08 0.16 0.35 0.08]);
plot(hAxX, 1:num_pixels, illum_gauss(num_pixels/2,:), 'k');
set(hAxX, 'Box','off', 'XLim', [0 1000]);

linkaxes([hAxImg hAxX], 'x');
linkaxes([hAxImg hAxY], 'y');

% Intercept comparison plot
hAxInt = axes(f,'Position',[0.55 0.25 0.40 0.65]);
hold(hAxInt,'on');
hPlotG = plot(hAxInt, nan, nan, 'b-',  'LineWidth',1.4, 'DisplayName','Gaussian');
hPlotF = plot(hAxInt, nan, nan, 'r--', 'LineWidth',1.4, 'DisplayName','Flat-top');
grid(hAxInt,'on');
xlabel(hAxInt,'Iteration');
ylabel(hAxInt,'BRDF / P_{tot}');
title(hAxInt,'Blue = Gaussian    Red = Flat-top');
legend(hAxInt,'Location','northwest');

%% 6. MAIN LOOP --------------------------------------------------------------
yG = nan(max_iter,1,'single');
yF = nan(max_iter,1,'single');

for it = 1:max_iter
    % Grow & brighten defects in both maps
    for k = 1:num_scatterers
        sc_r(k) = sc_r(k) + growth_mm;
        mask = ((x_mm - sc_cx(k)).^2 + (y_mm - sc_cy(k)).^2) ...
               <= (sc_r(k) * px_per_mm)^2;
        scatterG(mask) = scatterG(mask) + bright_step;
        scatterF(mask) = scatterF(mask) + bright_step;
    end
    
    % Compute weighted images
    imgG = illum_gauss .* scatterG + dark_level * (disk_mask & ~beam_mask);
    imgF = illum_flat  .* scatterF + dark_level * (disk_mask & ~beam_mask);
    
    % Update live scatter display (Gaussian branch)
    set(hImg,'CData',scatterG);
    
    % Compute ROI intercepts
    sumsG = zeros(numel(ROI_radii_mm),1,'single');
    sumsF = zeros(size(sumsG),'single');
    areas = zeros(size(sumsG),'single');
    for i = 1:numel(ROI_radii_mm)
        m          = ROI_masks(:,:,i);
        sumsG(i)   = sum(imgG(m));
        sumsF(i)   = sum(imgF(m));
        areas(i)   = sum(m(:));
    end
    
    % Linear fits → intercepts → normalisation
    coeffG    = polyfit(double(areas), double(sumsG), 1);
    coeffF    = polyfit(double(areas), double(sumsF), 1);
    yG(it)    = coeffG(2) / Ptot_gauss;
    yF(it)    = coeffF(2) / Ptot_flat;
    
    % Update plots
     set(hPlotG, 'XData', 1:it, 'YData', yG(1:it));
     set(hPlotF, 'XData', 1:it, 'YData', yF(1:it));
     legend(hAxInt, 'off');
    
    drawnow;
end