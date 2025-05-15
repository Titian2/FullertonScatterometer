%% Plot 1: All y-intercept curves by experiment class
%----------------------------------------------------

% 1) Read summary to get grouping & isNegative stats
opts = detectImportOptions('allExperiments_summary.tsv', ...
    'FileType','text','Delimiter','\t');
T = readtable('allExperiments_summary.tsv', opts);
% extract experiment class prefix before first underscore
grp = cellfun(@(s) strtok(s,'_'), T.testName, 'UniformOutput',false);
T.group = categorical(grp);
classes = categories(T.group);

% 2) Make 2×2 subplots
figure('Name','Y-Intercept Curves by Class','NumberTitle','off');
for ci = 1:numel(classes)
    cls = classes{ci};
    subplot(2,2,ci);
    hold on;
    
    % load that class's curves
    fileName = sprintf('allExperiments_%sCurves.mat', lower(cls));
    data = load(fileName,'allCurves');
    curves = data.allCurves;
    
    % plot each run with color based on isNegative
    nRuns = numel(curves);
    idx = (T.group == cls); % indices for this class
    isNeg = T.isNegative(idx); % isNegative values for this class
    
    for k = 1:nRuns
        if isNeg(k)
            plot(curves{k}, 'Color', 'r', 'LineWidth', 1); % red for negative
        else
            plot(curves{k}, 'Color', 'b', 'LineWidth', 1); % blue for positive
        end
    end
    
    % add legend for the first plot
    if ci == 1
        legend({'Negative', 'Positive'}, 'Location', 'best');
    end
     
    % compute how many runs ended negative
    idx = (T.group == cls);
    nNeg = sum(T.isNegative(idx));
    title(sprintf('%s  (neg: %d/%d)', cls, nNeg, sum(idx)), 'Interpreter','none');
    xlabel('Iteration'); ylabel('Y-Intercept');
    grid on; hold off;
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