%% Plot 1: All y-intercept curves by experiment class
%----------------------------------------------------


% Update file paths to include the new directory
basePath = '/Users/simon/Library/CloudStorage/GoogleDrive-simon.tait@ligo.org/My Drive/BackupFromDropbox/Fullerton_AAS/FullertonScatterometer/FullertonScatterometer/ModelResults/';
opts = detectImportOptions(fullfile(basePath, 'allExperiments_summary.tsv'), ...
    'FileType','text','Delimiter','\t');
T = readtable(fullfile(basePath, 'allExperiments_summary.tsv'), opts);

% extract experiment class prefix before first underscore
grp = cellfun(@(s) strtok(s,'_'), T.testName, 'UniformOutput',false);
T.group = categorical(grp);
classes = categories(T.group);

% === Plot 1: All curves by class ===
figure('Name','Y-Intercept Curves by Class','NumberTitle','off');

for ci = 1:numel(classes)
    cls = classes{ci};
    idxRows = find(T.group == cls);         % summary rows for this class
    nRows   = numel(idxRows);
    
    % Load curves for this class
    matFile = sprintf('allExperiments_%sCurves.mat', lower(cls));
    data = load(matFile, 'allCurves');
    curves = data.allCurves;                 % cell array
    nCurves = numel(curves);
    
    % How many to plot?
    nPlot = min(nRows, nCurves);
    if nPlot < nRows || nPlot < nCurves
      warning('Class %s: %d rows vs %d curves → plotting %d', ...
        cls, nRows, nCurves, nPlot);
    end
    
    % Subplot
    ax = subplot(2,2,ci); hold(ax,'on');    
    cmap = lines(nPlot);
    for k = 1:nPlot
      plot(ax, curves{k}, 'Color', cmap(k,:), 'LineWidth', 1);
    end
    
    % Negative count
    negCount = sum(T.isNegative(idxRows));
    title(ax, sprintf('%s   (neg: %d/%d)', cls, negCount, nRows), ...
          'Interpreter','none');
    xlabel(ax,'Iteration'); ylabel(ax,'Y-Intercept');
    grid(ax,'on');
end


%% Plot 2: Final-intercept vs. class, colored by isNegative
%---------------------------------------------------------

% numeric x‐positions for each class
cats = categories(T.group);
xpos = zeros(height(T),1);
for i = 1:numel(cats)
    xpos(T.group==cats{i}) = i;
end

% jitter to separate points
rng(0); 
j = (rand(size(xpos))-.5)*0.2;
xj = xpos + j;

% scatter non-negative (blue) and negative (red)
figure('Name','Final Intercept by Class & Outcome','NumberTitle','off');
hold on;
mask0 = T.isNegative==0;
mask1 = T.isNegative==1;
scatter(xj(mask0), T.finalInt(mask0), 20, 'b', 'filled');
scatter(xj(mask1), T.finalInt(mask1), 20, 'r', 'filled');
hold off;

xlim([0.5, numel(cats)+0.5]);
xticks(1:numel(cats));
xticklabels(cats);
xtickangle(45);
ylabel('Final Y-Intercept');
xlabel('Experiment Class');
title('Final Intercept by Class & isNegative');
legend({'Non-Negative','Negative'}, 'Location','best');
grid on;