%% Model for why the TIS estimation works — Circular ROI version
% Adapted from Csuflog20140820 toy model

clear all;
close all;
clc;

% 1) Build the toy image
background = ones(100,100);
embed      = ones(10,10)*10 + 1;
[l_1, l_2] = size(embed);
% center embed at (50,50)
% y_start = 50 - floor(l_1/2);
% x_start = 50 - floor(l_2/2);
% background(y_start:y_start+l_1-1, x_start:x_start+l_2-1) = embed;
img = background;

% 2) Define center of all ROIs
xcenter = 50;
ycenter = 50;

% 3) Base width → used to set circle radii
width_s = l_1 * 1.1;         % ~11 pixels
radii   = width_s * [0.5, 1, 1.5, 2];  % radii for circles

% 4) Precompute coordinate grid
[M, N] = size(img);
[X, Y] = meshgrid(1:N, 1:M);

% 5) Preallocate arrays
numROIs = numel(radii);
sums    = zeros(1, numROIs);
areas   = zeros(1, numROIs);

% 6) Loop over each circle ROI
for i = 1:numROIs
    mask      = (X - xcenter).^2 + (Y - ycenter).^2 <= radii(i).^2;
    sums(i)   = sum( img(mask) );
    areas(i)  = sum( mask(:) );    % <-- correct parentheses here
end

% 7) Fit line to sum vs. area
coeffs  = polyfit(areas, sums, 1);
m       = coeffs(1);
b       = coeffs(2);
fitLine = m * areas + b;

% 8) Show image with circles overlaid
figure;
imagesc(img);
colormap('gray');
axis image off;
hold on;
theta  = linspace(0, 2*pi, 360);
colors = {'b--', 'r--', 'g--', 'y--'};
for i = 1:numROIs
    xC = xcenter + radii(i) * cos(theta);
    yC = ycenter + radii(i) * sin(theta);
    plot(xC, yC, colors{i}, 'LineWidth', 2);
end
title('Circular ROIs overlaid','FontSize',14);

% 9) Plot sums vs. area and display intercept
figure;
plot(areas, sums, 'k^','MarkerSize',8,'LineWidth',1.5);
hold on;
plot(areas, fitLine, 'r-','LineWidth',2);
xlabel('ROI area (pixels)','FontSize',14);
ylabel('Sum of counts','FontSize',14);
title(sprintf('Sum vs. Area, y-intercept = %.1f counts', b),'FontSize',16);
legend('Measured','Linear fit','Location','NorthWest');
grid on;


%% 
%% Model for why the TIS estimation works — Circular ROIs + Random-circle background
% Adapted from Csuflog20140820 toy model

clear; close all; clc;

% 1) Build the toy image: background of random circles + central embed
imgSize = 100;
background = ones(imgSize, imgSize);

% grid for drawing circles
[Xbg, Ybg] = meshgrid(1:imgSize, 1:imgSize);

% scatter a bunch of random circles
numRandomCircles = 20;
minRadius = 3; 
maxRadius = 15;
for k = 1:numRandomCircles
    cx = randi([1, imgSize]);           % random center x
    cy = randi([1, imgSize]);           % random center y
    r  = randi([minRadius, maxRadius]); % random radius
    mask = (Xbg - cx).^2 + (Ybg - cy).^2 <= r.^2;
    % add a random extra count (0–10) inside each circle
    background(mask) = background(mask) + rand()*10;
end

% now embed the 10×10 high-scatter “spot” in the center
embed = ones(10,10)*10 + 1;
[l1, l2] = size(embed);
yStart = round(imgSize/2 - l1/2);
xStart = round(imgSize/2 - l2/2);
background(yStart:yStart+l1-1, xStart:xStart+l2-1) = embed;

img = background;


% 2) Define center for circular ROIs
xcenter = imgSize/2;
ycenter = imgSize/2;

% 3) Base width → determine ROI radii
width_s = l1 * 1.1;  % ~11 pixels
radii   = width_s * [0.5, 1, 1.5, 2];

% 4) Build coordinate grid for ROI masks
[M, N] = size(img);
[X, Y] = meshgrid(1:N, 1:M);

% 5) Preallocate sums & areas
numROIs = numel(radii);
sums    = zeros(1, numROIs);
areas   = zeros(1, numROIs);

% 6) Compute total counts & pixel-counts within each circular ROI
for i = 1:numROIs
    mask      = (X - xcenter).^2 + (Y - ycenter).^2 <= radii(i).^2;
    sums(i)   = sum( img(mask) );
    areas(i)  = sum( mask(:) );    % <-- fixed parentheses here
end

% 7) Fit sum vs. area to a line: sum = m·area + b
coeffs  = polyfit(areas, sums, 1);
m       = coeffs(1);
b       = coeffs(2);
fitLine = m * areas + b;

% 8) Show image with circular ROIs overlaid
figure;
imagesc(img);
colormap('gray');
axis image off;
hold on;
theta  = linspace(0, 2*pi, 360);
colors = {'b--','r--','g--','y--'};
for i = 1:numROIs
    xC = xcenter + radii(i)*cos(theta);
    yC = ycenter + radii(i)*sin(theta);
    plot(xC, yC, colors{i}, 'LineWidth', 2);
end
title('Circular ROIs overlaid','FontSize',14);

% 9) Plot sums vs. area and annotate intercept
figure;
plot(areas, sums, 'k^','MarkerSize',8,'LineWidth',1.5);
hold on;
plot(areas, fitLine, 'r-','LineWidth',2);
xlabel('ROI area (pixels)','FontSize',14);
ylabel('Sum of counts','FontSize',14);
title(sprintf('Sum vs. Area, y-intercept = %.1f counts', b),'FontSize',16);
legend('Measured','Linear fit','Location','NorthWest');
grid on;

%%

%% Model for why the TIS estimation works — 6 Nested Circular ROIs
% Adapted from Csuflog20140820 toy model

clear; close all; clc;

% 1) Build the toy image: background of random circles + central embed
imgSize = 100;
background = ones(imgSize, imgSize);

[Xbg, Ybg] = meshgrid(1:imgSize, 1:imgSize);
numRandomCircles = 20;
minRadius = 3; maxRadius = 15;
for k = 1:numRandomCircles
    cx = randi(imgSize);
    cy = randi(imgSize);
    r  = randi([minRadius, maxRadius]);
    mask = (Xbg - cx).^2 + (Ybg - cy).^2 <= r.^2;
    background(mask) = background(mask) + rand()*10;
end

embed = ones(10,10)*10 + 1;
[l1, l2] = size(embed);
yStart = round(imgSize/2 - l1/2);
xStart = round(imgSize/2 - l2/2);
background(yStart:yStart+l1-1, xStart:xStart+l2-1) = embed;

img = background;

% 2) Center for ROI placement
xcenter = imgSize/2;
ycenter = imgSize/2;

% 3) Define six nested radii, each 0.5*width_s apart
width_s = l1 * 1.1;           % ~11 pixels
factors = linspace(2.5,3,6);        % [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
radii   = width_s * factors;  % six radii

% 4) Precompute grid for mask creation
[M, N] = size(img);
[X, Y] = meshgrid(1:N, 1:M);

% 5) Allocate storage
numROIs = numel(radii);
sums    = zeros(1, numROIs);
areas   = zeros(1, numROIs);

% 6) Compute sums & areas inside each circular ROI
for i = 1:numROIs
    mask     = (X - xcenter).^2 + (Y - ycenter).^2 <= radii(i).^2;
    sums(i)  = sum( img(mask) );
    areas(i) = sum( mask(:) );
end

% 7) Linear fit: sum ≈ m·area + b
coeffs  = polyfit(areas, sums, 1);
m       = coeffs(1);
b       = coeffs(2);
fitLine = m * areas + b;

% 8) Display image with the six circles overlaid
figure;
imagesc(img);
colormap('gray');
axis image off; hold on;
theta  = linspace(0,2*pi,360);
colors = {'b--','r--','g--','y--','m--','c--'};
for i = 1:numROIs
    xC = xcenter + radii(i)*cos(theta);
    yC = ycenter + radii(i)*sin(theta);
    plot(xC, yC, colors{i}, 'LineWidth', 2);
end
title('Six Nested Circular ROIs','FontSize',14);

% 9) Plot sum vs. area with fit
figure;
plot(areas, sums, 'k^','MarkerSize',8,'LineWidth',1.5);
hold on;
plot(areas, fitLine, 'r-','LineWidth',2);
xlabel('ROI area (pixels)','FontSize',14);
ylabel('Sum of counts','FontSize',14);
title(sprintf('Sum vs. Area, y-intercept = %.1f counts', b),'FontSize',16);
legend('Measured','Linear fit','Location','NorthWest');
grid on;

%% 

%% Dynamic TIS intercept vs. background brightness with live subplot
% 6 nested circular ROIs, ramp the background over 1000 iterations,
% update both the ROI image and a live plot of the y-intercept evolution.

clear; close all; clc;

% 1) Build the toy image: random circles + central embed
imgSize = 100;
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

% embed = ones(10,10)*10 + 1;
% [l1, ~] = size(embed);
% yStart = round(imgSize/2 - l1/2);
% xStart = round(imgSize/2 - l1/2);
% background(yStart:yStart+l1-1, xStart:xStart+l1-1) = embed;

img = background;

% 2) Define 6 nested circular ROIs
l1 = 10 ;
xcenter = imgSize/2;
ycenter = imgSize/2;
width_s = l1 * 1.1;
radii   = width_s * 2.5 % linspace(2.5,3,6);  % six radii
[X, Y] = meshgrid(1:imgSize, 1:imgSize);
numROIs = numel(radii);

% 3) Preallocate storage for fit intercepts
numIterations = 1000;
yIntercepts   = nan(1, numIterations);

% 4) Create figure with two subplots
hFig = figure('Name','Dynamic ROIs and Intercept','NumberTitle','off',...
              'Position',[200 200 900 400]);

% Left: image + ROIs
hAxImg = subplot(1,2,1);
hImg = imagesc(img, 'Parent', hAxImg);
colormap(hAxImg, 'gray');
axis(hAxImg, 'image', 'off');
hold(hAxImg, 'on');
theta = linspace(0,2*pi,360);
colors = {'b--','r--','g--','y--','m--','c--'};
for i = 1:numROIs
    xC = xcenter + radii(i)*cos(theta);
    yC = ycenter + radii(i)*sin(theta);
    plot(hAxImg, xC, yC, colors{i}, 'LineWidth', 2);
end
title(hAxImg, 'Nested Circular ROIs','FontSize',14);

% Right: live plot of y-intercept vs iteration
hAxPlot = subplot(1,2,2);
hPlot = plot(hAxPlot, nan, nan, 'b-','LineWidth',1.5);
xlabel(hAxPlot, 'Iteration','FontSize',12);
ylabel(hAxPlot, 'Y-Intercept (counts)','FontSize',12);
title(hAxPlot, 'Evolution of Y-Intercept','FontSize',14);
xlim(hAxPlot, [1 numIterations]);
grid(hAxPlot, 'on');
hold(hAxPlot, 'on');

% 5) Main loop: ramp brightness, recalc intercept, update both plots
for iter = 1:numIterations
    % Increase brightness inside random circles
    img(randomMask) = img(randomMask) + 1;
    
    % Compute sums & areas for each ROI
    sums  = zeros(1, numROIs);
    areas = zeros(1, numROIs);
    for i = 1:numROIs
        maskROI   = (X - xcenter).^2 + (Y - ycenter).^2 <= radii(i).^2;
        sums(i)   = sum( img(maskROI) );
        areas(i)  = sum( maskROI(:) );
    end
    
    % Fit line: sum ≈ m·area + b
    coeffs           = polyfit(areas, sums, 1);
    yIntercepts(iter)= coeffs(2);
    
    % Update image data
    set(hImg, 'CData', img);
    
    % Update live intercept plot
    set(hPlot, 'XData', 1:iter, ...
               'YData', yIntercepts(1:iter));
    
    drawnow;
end

%%

%% Dynamic TIS intercept vs. growing & brightening circles
% 6 nested circular ROIs, with random background circles that start small
% (5 px) and then grow in size and brightness over 1000 iterations; live
% update of ROI image + intercept plot.

clear; close all; clc;

% 1) Initialize toy image with random circles0
imgSize = 4096;
background = ones(imgSize);  % base level = 1

% Precompute grid for drawing circles
[Xbg, Ybg] = meshgrid(1:imgSize, 1:imgSize);

% Generate random circle parameters
numRandomCircles = 100;
initialRadius    = 123;                        % start every circle at 5 px
cx    = randi(imgSize, [1, numRandomCircles]);  
cy    = randi(imgSize, [1, numRandomCircles]);
r0    = ones(1, numRandomCircles) * initialRadius;  
initAmp = rand(1, numRandomCircles) * 10;    % initial extra amplitude

% Paint initial random circles (all radius = 5 px)
for k = 1:numRandomCircles
    mask0 = (Xbg - cx(k)).^2 + (Ybg - cy(k)).^2 <= r0(k)^2;
    background(mask0) = background(mask0) + initAmp(k);
end

img = background;  % working image

% 2) Define 6 nested circular ROIs (fixed)
xcenter = imgSize/2;
ycenter = imgSize/2;
width_s = 400 * 1.1;                      
radii   = width_s * linspace(2.5, 3, 6); % six radii, closely spaced

% Precompute grid for ROI masks
[X, Y] = meshgrid(1:imgSize, 1:imgSize);
numROIs = numel(radii);

% 3) Preallocate storage for y-intercepts
numIterations = 1000;
yIntercepts   = nan(1, numIterations);

% 4) Set growth parameters for random circles
growthRate     = width_s*0.01;    % px per iteration
brightnessStep = 0.1;     % counts added per iteration inside circles

% 5) Create figure with two subplots
hFig = figure('Name','Growing Circles & ROI Fit','NumberTitle','off',...
              'Position',[200 200 900 400]);

% Left subplot: ROI image
hAxImg = subplot(1,2,1);
hImg   = imagesc(img, 'Parent', hAxImg);
colormap(hAxImg,'gray');
axis(hAxImg,'image','off');
hold(hAxImg,'on');
theta    = linspace(0,2*pi,360);
colorsROI = {'b--','r--','g--','y--','m--','c--'};
for i = 1:numROIs
    xC = xcenter + radii(i)*cos(theta);
    yC = ycenter + radii(i)*sin(theta);
    plot(hAxImg, xC, yC, colorsROI{i}, 'LineWidth', 2);
end
title(hAxImg,'Nested Circular ROIs','FontSize',14);

% Right subplot: live y-intercept plot
hAxPlot = subplot(1,2,2);
hPlot   = plot(hAxPlot, nan, nan, 'b-','LineWidth',1.5);
xlabel(hAxPlot,'Iteration','FontSize',12);
ylabel(hAxPlot,'Y-Intercept (counts)','FontSize',12);
title(hAxPlot,'Evolution of Y-Intercept','FontSize',14);
xlim(hAxPlot,[1 numIterations]);
grid(hAxPlot,'on');
hold(hAxPlot,'on');

% 6) Main loop: grow & brighten circles, recalc fit, update plots
for iter = 1:numIterations
    % Build dynamic mask for all circles at their current radii
    currentMask = false(imgSize);
    for k = 1:numRandomCircles
        rCurr    = r0(k) + growthRate * iter;
        thisMask = (Xbg - cx(k)).^2 + (Ybg - cy(k)).^2 <= rCurr^2;
        currentMask = currentMask | thisMask;
    end
    
    % Increase brightness inside growing circles
    img(currentMask) = img(currentMask) + brightnessStep;
    
    % Compute ROI sums & areas
    sums  = zeros(1, numROIs);
    areas = zeros(1, numROIs);
    for i = 1:numROIs
        maskROI   = (X - xcenter).^2 + (Y - ycenter).^2 <= radii(i)^2;
        sums(i)   = sum(img(maskROI));
        areas(i)  = sum(maskROI(:));  % correct parentheses
    end
    
    % Linear fit: sums ≈ m·areas + b
    coeffs            = polyfit(areas, sums, 1);
    yIntercepts(iter) = coeffs(2);
    
    % Update ROI image
    set(hImg, 'CData', img);
    
    % Update live intercept plot
    set(hPlot, 'XData', 1:iter, 'YData', yIntercepts(1:iter));
    
    drawnow;
end

% After loop, yIntercepts contains the evolution of the background intercept.



%% Actual Image Size - Slow 

clear; close all; clc;

% 1) Initialize toy image with random circles
imgSize = 4096;
background = ones(imgSize)*1E-5;

[Xbg, Ybg] = meshgrid(1:imgSize, 1:imgSize);
numRandomCircles = 100;
initialRadius    = 123;
cx = randi(imgSize, [1,numRandomCircles]);
cy = randi(imgSize, [1,numRandomCircles]);
r0 = ones(1,numRandomCircles)*initialRadius;
initAmp = rand(1,numRandomCircles)*10;

for k = 1:numRandomCircles
    mask0 = (Xbg-cx(k)).^2 + (Ybg-cy(k)).^2 <= r0(k)^2;
    background(mask0) = background(mask0) + initAmp(k);
end

img = background;

% 2) Define nested circular ROIs
xcenter = imgSize/2;
ycenter = imgSize/2;
width_s = 400 * 1.1;
radii   = width_s * linspace(2.5,3,6);
[X, Y]  = meshgrid(1:imgSize,1:imgSize);
numROIs = numel(radii);

% 3) Preallocate
numIterations = 1000;
yIntercepts   = nan(1,numIterations);

% 4) Growth params
growthRate     = width_s*0.01;
brightnessStep = 0.1;

% 5) Create figure with 3 subplots
hFig = figure('Name','Growing Circles & ROI Fit','NumberTitle','off',...
              'Position',[200 200 1200 400]);

% a) ROI image
hAxImg = subplot(1,3,1);
hImg = imagesc(img,'Parent',hAxImg);
colormap(hAxImg,'gray');
axis(hAxImg,'image','off'); hold(hAxImg,'on');
theta = linspace(0,2*pi,360);
colorsROI = {'b--','r--','g--','y--','m--','c--'};
for i=1:numROIs
    xC = xcenter + radii(i)*cos(theta);
    yC = ycenter + radii(i)*sin(theta);
    plot(hAxImg, xC, yC, colorsROI{i}, 'LineWidth',2);
end
title(hAxImg,'Nested Circular ROIs','FontSize',14);

% b) Y-intercept vs iteration
hAxPlot = subplot(1,3,2);
hPlot   = plot(hAxPlot, nan, nan, 'b-','LineWidth',1.5);
xlabel(hAxPlot,'Iteration','FontSize',12);
ylabel(hAxPlot,'Y-Intercept','FontSize',12);
title(hAxPlot,'Evolution of Y-Intercept','FontSize',14);
xlim(hAxPlot,[1 numIterations]);
grid(hAxPlot,'on'); hold(hAxPlot,'on');

% c) Sums vs area + fit
hAxSA   = subplot(1,3,3);
hSA_pts = plot(hAxSA, nan, nan, 'k^','MarkerSize',8,'LineWidth',1.5);
hold(hAxSA,'on');
hSA_line= plot(hAxSA, nan, nan, 'r-','LineWidth',2);
xlabel(hAxSA,'ROI area (pixels)','FontSize',14);
ylabel(hAxSA,'Sum of counts','FontSize',14);
title(hAxSA,'Sum vs. Area','FontSize',14);
legend(hAxSA,'Measured','Linear fit','Location','NorthWest');
grid(hAxSA,'on');

% 6) Main loop
for iter = 1:numIterations
    % update circles
    currentMask = false(imgSize);
    for k = 1:numRandomCircles
        rCurr = r0(k) + growthRate*iter;
        m = (Xbg-cx(k)).^2 + (Ybg-cy(k)).^2 <= rCurr^2;
        currentMask = currentMask | m;
    end
    img(currentMask) = img(currentMask) + brightnessStep;

    % compute sums & areas
    sums  = zeros(1,numROIs);
    areas = sums;
    for i = 1:numROIs
        maskROI   = (X-xcenter).^2 + (Y-ycenter).^2 <= radii(i)^2;
        sums(i)   = sum(img(maskROI));
        areas(i)  = sum(maskROI(:));
    end

    % fit line
    coeffs            = polyfit(areas, sums, 1);
    m                 = coeffs(1);
    b                 = coeffs(2);
    fitLine           = m*areas + b;
    yIntercepts(iter) = b;

    % refresh plots
    set(hImg,   'CData', img);
    set(hPlot, 'XData', 1:iter,       'YData', yIntercepts(1:iter));
    set(hSA_pts,'XData', areas,       'YData', sums);
    set(hSA_line,'XData', areas,      'YData', fitLine);
    title(hAxSA, sprintf('Sum vs. Area (b=%.1f)', b), 'FontSize',14);

    drawnow;
end


%% USING Paralell processing

clear; close all; clc;

% 1) Initialize toy image with random circles
imgSize = 4096;
background = ones(imgSize);

[Xbg, Ybg] = meshgrid(1:imgSize, 1:imgSize);
numRandomCircles = 100;
initialRadius    = 123;
cx = randi(imgSize, [1,numRandomCircles]);
cy = randi(imgSize, [1,numRandomCircles]);
r0 = ones(1,numRandomCircles)*initialRadius;
initAmp = rand(1,numRandomCircles)*10;

for k = 1:numRandomCircles
    mask0 = (Xbg-cx(k)).^2 + (Ybg-cy(k)).^2 <= r0(k)^2;
    background(mask0) = background(mask0) + initAmp(k);
end

img = background;

% 2) Define nested circular ROIs
xcenter = imgSize/2;
ycenter = imgSize/2;
width_s = 400 * 1.1;
radii   = width_s * linspace(2.5,3,6);
[X, Y]  = meshgrid(1:imgSize,1:imgSize);
numROIs = numel(radii);

% 3) Preallocate
numIterations = 1000;
yIntercepts   = nan(1,numIterations);

% 4) Growth params
growthRate     = width_s*0.01;
brightnessStep = 0.1;

% 5) Prepare figure with 3 subplots
hFig = figure('Name','Growing Circles & ROI Fit','NumberTitle','off',...
              'Position',[200 200 1200 400]);

% a) ROI image
hAxImg = subplot(1,3,1);
hImg = imagesc(img,'Parent',hAxImg);
colormap(hAxImg,'gray');
axis(hAxImg,'image','off'); hold(hAxImg,'on');
theta = linspace(0,2*pi,360);
colorsROI = {'b--','r--','g--','y--','m--','c--'};
for i = 1:numROIs
    xC = xcenter + radii(i)*cos(theta);
    yC = ycenter + radii(i)*sin(theta);
    plot(hAxImg, xC, yC, colorsROI{i}, 'LineWidth',2);
end
title(hAxImg,'Nested Circular ROIs','FontSize',14);

% b) Y-intercept vs iteration
hAxPlot = subplot(1,3,2);
hPlot   = plot(hAxPlot, nan, nan, 'b-','LineWidth',1.5);
xlabel(hAxPlot,'Iteration','FontSize',12);
ylabel(hAxPlot,'Y-Intercept','FontSize',12);
title(hAxPlot,'Evolution of Y-Intercept','FontSize',14);
xlim(hAxPlot,[1 numIterations]);
grid(hAxPlot,'on'); hold(hAxPlot,'on');

% c) Sums vs area + fit
hAxSA    = subplot(1,3,3);
hSA_pts  = plot(hAxSA, nan, nan, 'k^','MarkerSize',8,'LineWidth',1.5);
hold(hAxSA,'on');
hSA_line = plot(hAxSA, nan, nan, 'r-','LineWidth',2);
xlabel(hAxSA,'ROI area (pixels)','FontSize',14);
ylabel(hAxSA,'Sum of counts','FontSize',14);
title(hAxSA,'Sum vs. Area','FontSize',14);
legend(hAxSA,'Measured','Linear fit','Location','NorthWest');
grid(hAxSA,'on');

% 6) Main loop with parfor mask build
for iter = 1:numIterations
    % Build mask indices in parallel
    idxList = cell(numRandomCircles,1);
    parfor k = 1:numRandomCircles
        rCurr = r0(k) + growthRate * iter;
        localMask = (Xbg-cx(k)).^2 + (Ybg-cy(k)).^2 <= rCurr^2;
        idxList{k} = find(localMask);
    end
    % Union all indices into one mask
    allIdx = vertcat(idxList{:});
    maskAll = false(imgSize);
    maskAll(allIdx) = true;

    % Increase brightness inside growing circles
    img(maskAll) = img(maskAll) + brightnessStep;

    % Compute ROI sums & areas
    sums  = zeros(1,numROIs);
    areas = sums;
    for i = 1:numROIs
        maskROI   = (X-xcenter).^2 + (Y-ycenter).^2 <= radii(i)^2;
        sums(i)   = sum(img(maskROI));
        areas(i)  = sum(maskROI(:));
    end

    % Fit line
    coeffs            = polyfit(areas, sums, 1);
    m                 = coeffs(1);
    b                 = coeffs(2);
    fitLine           = m*areas + b;
    yIntercepts(iter) = b;

    % Refresh plots
    set(hImg,    'CData', img);
    set(hPlot,  'XData', 1:iter,        'YData', yIntercepts(1:iter));
    set(hSA_pts,'XData', areas,        'YData', sums);
    set(hSA_line,'XData', areas,       'YData', fitLine);
    title(hAxSA, sprintf('Sum vs. Area (b=%.1f)', b), 'FontSize',14);

    drawnow;
end



