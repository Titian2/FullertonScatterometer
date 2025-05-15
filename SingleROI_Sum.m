%% Single ROI With SUM instead of Y_intercept 

clear; close all; clc;

% 1) Build the toy image
imgSize    = 100;
background = ones(imgSize);

[Xbg, Ybg] = meshgrid(1:imgSize, 1:imgSize);
numRandomCircles = 20;
minRadius = 3; maxRadius = 15;
randomMask = false(imgSize);

for k = 1:numRandomCircles
    cx = randi(imgSize);
    cy = randi(imgSize);
    r  = randi([minRadius, maxRadius]);
    thisMask = (Xbg - cx).^2 + (Ybg - cy).^2 <= r.^2;
    randomMask = randomMask | thisMask;
    background(thisMask) = background(thisMask) + rand()*10;
end

img = background;  % working image

% 2) Define a single circular ROI
l1      = 10;
xcenter = imgSize/2;
ycenter = imgSize/2;
width_s = l1 * 1.1;
radius  = width_s * 2.5;         % just one ROI now
[X, Y]   = meshgrid(1:imgSize, 1:imgSize);

% precompute ROI mask
maskROI = (X - xcenter).^2 + (Y - ycenter).^2 <= radius^2;

% 3) Preallocate storage for sum-in-ROI
numIterations = 1000;
sumROI        = zeros(1, numIterations);

% 4) Growth parameters
growthRate     = 1;      % your chosen value
brightnessStep = 1;      % your chosen value

% 5) Create figure with two subplots
hFig = figure('Name','Dynamic ROI Sum','NumberTitle','off',...
              'Position',[200 200 900 400]);

% a) Left: image
hAxImg = subplot(1,2,1);
hImg   = imagesc(img, 'Parent', hAxImg);
colormap(hAxImg,'gray');
axis(hAxImg,'image','off');
hold(hAxImg,'on');
theta = linspace(0,2*pi,360);
plot(hAxImg, xcenter + radius*cos(theta), ...
             ycenter + radius*sin(theta), 'r--','LineWidth',2);
title(hAxImg,'Image + Single ROI','FontSize',14);

% b) Right: live sum plot
hAxPlot = subplot(1,2,2);
hPlot   = plot(hAxPlot, nan, nan, 'b-','LineWidth',1.5);
xlabel(hAxPlot,'Iteration','FontSize',12);
ylabel(hAxPlot,'Sum inside ROI','FontSize',12);
title(hAxPlot,'ROI Sum vs Iteration','FontSize',14);
xlim(hAxPlot, [1 numIterations]);
grid(hAxPlot,'on');
hold(hAxPlot,'on');

% 6) Main loop: ramp brightness, update sum plot
for iter = 1:numIterations
    % brighten randomly masked region
    img(randomMask) = img(randomMask) + brightnessStep;

    % compute sum inside the single ROI
    sumROI(iter) = sum( img(maskROI) );

    % update image display
    set(hImg,   'CData', img);

    % update live sum plot
    set(hPlot, 'XData', 1:iter, ...
               'YData', sumROI(1:iter));

    drawnow;
end

%% 

clear; close all; clc;

% 1) Build the toy image
imgSize    = 100;
background = ones(imgSize);

% Precompute ROI geometry
l1      = 10;
xcenter = imgSize/2;
ycenter = imgSize/2;
width_s = l1 * 1.1;
radius  = width_s * 2.5;
[X, Y]   = meshgrid(1:imgSize, 1:imgSize);
maskROI  = (X - xcenter).^2 + (Y - ycenter).^2 <= radius^2;

% 2) Generate random circles, but *only inside* that ROI
numRandomCircles = 20;
minR = 3; maxR = 15;
randomMask = false(imgSize);

for k = 1:numRandomCircles
    % sample a random point inside the ROI circle
    r_samp    = radius * sqrt(rand());      % sqrt for uniform area
    theta_samp= 2*pi * rand();
    cx = round(xcenter + r_samp * cos(theta_samp));
    cy = round(ycenter + r_samp * sin(theta_samp));
    % ensure inside bounds
    cx = min(max(cx,1),imgSize);
    cy = min(max(cy,1),imgSize);
    
    % now paint a small circle of random radius & amplitude
    rr = randi([minR, maxR]);
    thisMask = (X-cx).^2 + (Y-cy).^2 <= rr^2;
    randomMask = randomMask | thisMask;
    background(thisMask) = background(thisMask) + rand()*10;
end

img = background;  % working image

% 3) Preallocate storage for sum‐in‐ROI
numIterations = 1000;
sumROI        = zeros(1, numIterations);

% 4) Growth parameters
growthRate     = 1;
brightnessStep = 1;

% 5) Create figure with two subplots
hFig = figure('Name','Dynamic ROI Sum','NumberTitle','off',...
              'Position',[200 200 900 400]);

% a) Left: image + ROI
hAxImg = subplot(1,2,1);
hImg   = imagesc(img, 'Parent', hAxImg);
colormap(hAxImg,'gray');
axis(hAxImg,'image','off');
hold(hAxImg,'on');
theta = linspace(0,2*pi,360);
plot(hAxImg, xcenter + radius*cos(theta), ...
             ycenter + radius*sin(theta), 'r--','LineWidth',2);
title(hAxImg,'Image with ROI','FontSize',14);

% b) Right: live sum plot
hAxPlot = subplot(1,2,2);
hPlot   = plot(hAxPlot, nan, nan, 'b-','LineWidth',1.5);
xlabel(hAxPlot,'Iteration','FontSize',12);
ylabel(hAxPlot,'Sum inside ROI','FontSize',12);
title(hAxPlot,'ROI Sum vs Iteration','FontSize',14);
xlim(hAxPlot, [1 numIterations]);
grid(hAxPlot,'on');
hold(hAxPlot,'on');

% 6) Main loop: ramp brightness, update sum plot
for iter = 1:numIterations
    % brighten only inside those initial random spots
    img(randomMask) = img(randomMask) + brightnessStep;

    % compute sum inside our single ROI
    sumROI(iter) = sum( img(maskROI) );

    % update image display
    set(hImg, 'CData', img);

    % update live sum plot
    set(hPlot, 'XData', 1:iter, ...
               'YData', sumROI(1:iter));

    drawnow;
end