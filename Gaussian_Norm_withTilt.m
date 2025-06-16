%% --------------------------------------------------------------------------
%  Scattering-map simulator – Gaussian vs. flat-top with laser instability
%  v 0.4 
% --------------------------------------------------------------------------
clear; close all; clc;

%% 1. USER PARAMETERS --------------------------------------------------------
sample_diam_mm    = 25.4;                   % disk diameter (mm)
beam_FWHM_mm      = [5 5];                  % Gaussian FWHM (x,y) in mm
beam_center_mm    = [0 0];                  % beam offset wrt disk centre
tilt_deg          = 10;                     % left→right 1/d^2 fall-off (deg)
num_pixels        = 1024;                   % simulation grid N×N

num_scatterers    = 1000;                   % number of scatterers
initial_diam_mm   = 0.005;                  % starting scatter diameter (mm)
growth_mm         = initial_diam_mm * 0.001;% radius growth per iter (mm)
bright_step       = 0.10;                   % brightness boost per iter
max_iter          = 1000;                   % total iterations

beam_threshold    = 0.10;                   % ≥10% of peak → inside beam
dark_level        = 0;                      % background outside beam
ROI_radii_mm      = linspace(2.5,3.5,6);    % ROI radii in mm

%% 2. GRID & STATIC MAPS -----------------------------------------------------
px_per_mm  = num_pixels / sample_diam_mm;
[Xpx,Ypx]  = meshgrid(1:num_pixels, 1:num_pixels);
x_mm       = (Xpx - (num_pixels+1)/2) / px_per_mm;
y_mm       = (Ypx - (num_pixels+1)/2) / px_per_mm;

disk_mask = (x_mm.^2 + y_mm.^2) <= (sample_diam_mm/2)^2;

% Tilt map
d0         = 1000; %% Distance from camera to sample in mm 
tilt_map   = (d0 ./ (d0 + x_mm * sind(tilt_deg))).^2 .* disk_mask;

% Gaussian beam map
sig_x      = beam_FWHM_mm(1)/(2*sqrt(2*log(2)));
sig_y      = beam_FWHM_mm(2)/(2*sqrt(2*log(2)));
beam_map   = exp(-((x_mm-beam_center_mm(1)).^2)/(2*sig_x^2) ...
                 -((y_mm-beam_center_mm(2)).^2)/(2*sig_y^2));
beam_map(~disk_mask) = 0;

% Illumination profiles
illum_gauss = beam_map .* tilt_map;                  % true Gaussian × tilt
beam_mask   = illum_gauss >= beam_threshold * max(illum_gauss(:));
illum_flat  = double(beam_mask);                     % uniform inside beam

Ptot_gauss  = sum(illum_gauss(:));
Ptot_flat   = sum(illum_flat(:));

% Laser instability: 1% sinusoidal fluctuation, 100-iteration period
laser_instability = 1 + 0.01 * sin(2*pi*(1:max_iter)/100);

%% 3. INITIAL DEFECT MAP (inside beam only) ---------------------------------
rng default
beam_lin    = find(beam_mask);
pick_inds   = beam_lin(randi(numel(beam_lin),1,num_scatterers));
[py, px]    = ind2sub([num_pixels,num_pixels], pick_inds);
sc_cx       = (px - (num_pixels+1)/2) / px_per_mm;
sc_cy       = (py - (num_pixels+1)/2) / px_per_mm;
sc_r        = initial_diam_mm * ones(1,num_scatterers);
sc_amp      = rand(1,num_scatterers) * 20;

scatterG    = zeros(num_pixels,'single');  % Gaussian branch
scatterF    = zeros(num_pixels,'single');  % flat-top branch
for k = 1:num_scatterers
    mask = ((x_mm-sc_cx(k)).^2 + (y_mm-sc_cy(k)).^2) <= (sc_r(k)*px_per_mm)^2;
    scatterG(mask) = scatterG(mask) + sc_amp(k);
    scatterF(mask) = scatterF(mask) + sc_amp(k);
end

%% 4. ROI MASKS --------------------------------------------------------------
ROI_masks = false(num_pixels, num_pixels, numel(ROI_radii_mm));
for i = 1:numel(ROI_radii_mm)
    ROI_masks(:,:,i) = ((x_mm-beam_center_mm(1)).^2 + ...
                        (y_mm-beam_center_mm(2)).^2) <= (ROI_radii_mm(i))^2;
end

%% 5. FIGURE LAYOUT ----------------------------------------------------------
f = figure('Color','w','Units','normalized','Position',[0.05 0.1 0.9 0.75]);

% Main image (Gaussian-weighted scatter)
hAxImg = axes(f,'Position',[0.08 0.25 0.35 0.65]);
hImg   = imagesc(scatterG,'Parent',hAxImg);
axis(hAxImg,'image','off'); colormap(hAxImg,'bone');
caxis(hAxImg,[0 max(scatterG(:))+eps]); hold(hAxImg,'on');

% ROI rings & disk rim
theta = linspace(0,2*pi,360);
for i = 1:numel(ROI_radii_mm)
    r_px = ROI_radii_mm(i)*px_per_mm;
    plot(hAxImg,(num_pixels+1)/2+r_px*cos(theta), ...
                 (num_pixels+1)/2+r_px*sin(theta),'m--','LineWidth',1.2);
end
r_px_disk = (sample_diam_mm/2)*px_per_mm;
plot(hAxImg,(num_pixels+1)/2+r_px_disk*cos(theta), ...
             (num_pixels+1)/2+r_px_disk*sin(theta),'w--','LineWidth',1.1);
title(hAxImg,'Defects (bone colormap)');

% Y-profile axis (manual position, 0–1000 px)
hAxY = axes(f,'Position',[0.017 0.24 0.05 0.66]);
plot(hAxY, illum_gauss(:,num_pixels/2), 1:num_pixels, 'k');
set(hAxY,'YDir','reverse','Box','off','YLim',[0 1000]);

% X-profile axis (manual position, 0–1000 px)
hAxX = axes(f,'Position',[0.08 0.16 0.35 0.08]);
plot(hAxX, 1:num_pixels, illum_gauss(num_pixels/2,:), 'k');
set(hAxX,'Box','off','XLim',[0 1000]);

linkaxes([hAxImg hAxX],'x'); linkaxes([hAxImg hAxY],'y');

% Intercept comparison plot
hAxInt = axes(f,'Position',[0.51,0.32,0.44,0.57]); hold(hAxInt,'on');
hPlotG = plot(hAxInt, nan, nan, 'bo-','LineWidth',1.4,'DisplayName','Gaussian');
hPlotF = plot(hAxInt, nan, nan, 'ro-','LineWidth',1.4,'DisplayName','Flat-top');
grid(hAxInt,'on');
xlabel(hAxInt,'Iteration');
ylabel(hAxInt,'Intercept / P_{tot}');
title(hAxInt,'Blue = Gaussian    Red = Flat-top');
legend(hAxInt,'Location','northwest');

% Power fluctuation subplot
hAxPow = axes(f,'Position',[0.51,0.12,0.44,0.12]); hold(hAxPow,'on');
hPlotP = plot(hAxPow, nan, nan, 'ko-','LineWidth',1.2);
set(hAxPow, ...
    'Box','off', ...
    'XLim',[1 max_iter], ...        % make sure both axes cover the same range
    'YLim',[min(laser_instability) max(laser_instability)]);
ylabel(hAxPow,'Power');
xlabel(hAxPow,'Iteration');
title(hAxPow,'Laser Instability');

% --- LINK X-AXES ---
linkaxes([hAxInt, hAxPow], 'x');



%% 6. MAIN LOOP --------------------------------------------------------------
yG = nan(max_iter,1,'single');
yF = nan(max_iter,1,'single');

% setup for unstable laser 
jitter   = zeros(max_iter,1,'double');      % to track scaleFactor
unstableMean = zeros(max_iter,1,'double');  % to track mean(unstable_laser_sim)

%% Main loop: propagate scaled noise, update images & plots\\

for it = 1:max_iter
    % --- Grow & brighten defects ---
    for k = 1:num_scatterers
        sc_r(k) = sc_r(k) + growth_mm;
        mask    = ((x_mm-sc_cx(k)).^2 + (y_mm-sc_cy(k)).^2) <= (sc_r(k)*px_per_mm)^2;
        scatterG(mask) = scatterG(mask) + bright_step;
        scatterF(mask) = scatterF(mask) + bright_step;
    end
    
    % --- Laser instability scaling + noise ---
    baseScale = laser_instability(it);
    scaleFactor = baseScale * (1 + 0.1*randn);            % add Gaussian noise
    noise_std   = 0.01 * max(scatterG(:));                % control noise level
    
    % create a per‐pixel “unstable laser” map:
        %   – first term: overall scale with 10% frame‐to‐frame jitter
        %   – second term: per‐pixel Gaussian noise
    unstable_laser_sim = baseScale*(1 + 0.1*randn) ...           % scalar scale jitter
                      + noise_std*randn(size(scatterG));    % pixel‐wise noise


    Noise_illum_gauss = unstable_laser_sim .* illum_gauss;
    Noise_illum_flat = unstable_laser_sim .* illum_flat;
    % take raw illumination maps and scale with the simulated laser power
    % adding on the random noise - this should be replaced with an actual
    % stream of laser power data 

    imgG = Noise_illum_gauss.* scatterG + dark_level * (disk_mask & ~beam_mask);   

    imgF = Noise_illum_flat  .* scatterF+ dark_level * (disk_mask & ~beam_mask);


  % store for plotting
    jitter(it)       = scaleFactor;
    unstableMean(it) = mean(unstable_laser_sim(:));     % should be ≈ scaleFactor

    % --- Update live image ---
    set(hImg,'CData',scatterG);

    % --- Compute ROI sums for intercept ---
    sumsG = zeros(numel(ROI_radii_mm),1,'single');
    sumsF = sumsG;
    areas = sumsG;
    for i = 1:numel(ROI_radii_mm)
        m = ROI_masks(:,:,i);
        sumsG(i) = sum(imgG(m));
        sumsF(i) = sum(imgF(m));
        areas(i) = sum(m(:));
    end
    
    % --- Linear fit & normalization ---
    coeffG = polyfit(double(areas), double(sumsG), 1);
    coeffF = polyfit(double(areas), double(sumsF), 1);
    yG(it) = coeffG(2) / Ptot_gauss;
    yF(it) = coeffF(2) / Ptot_flat;
    




    % --- Update intercept & power plots ---
    set(hPlotG, 'XData',1:it,      'YData',yG(1:it));
    set(hPlotF, 'XData',1:it,      'YData',yF(1:it));
    set(hPlotP, 'XData',1:it,      'YData',jitter(1:it));
    

    
    
    
    drawnow;
end