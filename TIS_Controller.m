%% === DOE Driver with per‐Sweep Waitbars ===

% Delete old combined summary
summaryFile = 'allExperiments_summary.tsv';
if isfile(summaryFile)
  delete(summaryFile);
end

% --- Base params (common) ---
baseParams = struct();

baseParams.imgSize          = 100;
baseParams.numRandomCircles = 20;
baseParams.initialRadius    = 5;
baseParams.embedSize        = 10;
baseParams.numIterations    = 1000;
baseParams.growthRate       = 0.01;
baseParams.brightnessStep   = 0.1;
baseParams.initAmpDist      = struct('type','uniform','min',0,'max',10);
baseParams.placementBias    = struct('type','uniform');
baseParams.radiiFactors     = linspace(2.5,3,6);

%% 1) Placement‐Bias Sweep with proper sigma check
biasConfigs = { struct('type','uniform'), ...
                struct('type','gaussian','sigma',10), ...
                struct('type','gaussian','sigma',20), ...
                struct('type','gaussian','sigma',30), ...
                struct('type','gaussian','sigma',40), ...
                struct('type','gaussian','sigma',50) };

nConfigs   = numel(biasConfigs);
runsPerCfg = 100;
totalRuns  = nConfigs * runsPerCfg;
h1 = waitbar(0, 'Placement bias sweep: 0%');

baseParams.testName         = 'GaussianTest';

% preallocate curve storage
allCurves = cell(totalRuns,1);

runCount = 0;
for c = 1:nConfigs
  cfg = biasConfigs{c};
  % extract sigma safely
  if isfield(cfg,'sigma')
    sigmaVal = cfg.sigma;
  else
    sigmaVal = 0;
  end

  for iter = 1:runsPerCfg
    runCount = runCount + 1;

    params = baseParams;
    params.placementBias = cfg;

    % build a reproducible seed:
    %   use an ASCII sum of the type string + sigmaVal + iter
    asciiSum = sum(double(cfg.type));  
    params.seed = uint32(asciiSum + sigmaVal + iter);

    result = runDynamicTIS(params, false);

    % store the full curve
    allCurves{runCount} = result.yIntercepts;

    waitbar(runCount/totalRuns, h1, ...
      sprintf('Placement bias: %d / %d (%.1f%%)', runCount, totalRuns, 100*runCount/totalRuns));
  end
end
close(h1);

% save all curves to disk
save('allExperiments_gaussiantestCurves.mat', 'allCurves');

%% 2) Amplitude‐Distribution Sweep

ampConfigs = { struct('type','uniform','min',0,'max',5), ...
               struct('type','uniform','min',0,'max',10), ...
               struct('type','uniform','min',0,'max',20), ...
               struct('type','normal','mu',5,'sigma',2), ...
               struct('type','normal','mu',5,'sigma',5), ...
               struct('type','normal','mu',10,'sigma',5) };
nConfigs   = numel(ampConfigs);
totalRuns  = nConfigs * runsPerCfg;
h2 = waitbar(0, 'Amplitude distribution sweep: 0%');

baseParams.testName         = 'Amplitude‐Distribution';

% preallocate curve storage
allCurves = cell(totalRuns,1);

runCount = 0;
for c = 1:nConfigs
  cfg = ampConfigs{c};
  for iter = 1:runsPerCfg
    runCount = runCount + 1;
    params = baseParams;
    params.initAmpDist = cfg;
    params.seed        = uint32(sum(double(string(cfg.type))) + sum(struct2array(cfg)) + iter);
    result = runDynamicTIS(params, false);

    % store the full curve
    allCurves{runCount} = result.yIntercepts;

    waitbar(runCount/totalRuns, h2, ...
      sprintf('Amplitude dist: %d / %d (%.1f%%)', runCount, totalRuns, 100*runCount/totalRuns));
  end
end
close(h2);

% save all curves to disk
save('allExperiments_amplitude‐distributionCurves.mat', 'allCurves');

%% 3) Growth‐vs‐Brightness Sweep
growthRates  = [0.001, 0.005, 0.01, 0.02];
brightSteps  = [0.05, 0.1, 0.2, 0.5];
nConfigs     = numel(growthRates)*numel(brightSteps);
totalRuns    = nConfigs * runsPerCfg;
h3 = waitbar(0, 'Growth vs brightness sweep: 0%');

baseParams.testName         = 'Growth‐vs‐Brightness';
% preallocate curve storage
allCurves = cell(totalRuns,1);

runCount = 0;
for gr = growthRates
  for br = brightSteps
    for iter = 1:runsPerCfg
      runCount = runCount + 1;
      params = baseParams;
      params.growthRate     = gr;
      params.brightnessStep = br;
      params.seed           = uint32(gr*1e5 + br*1e2 + iter);
      result= runDynamicTIS(params, false);
       % store the full curve
      allCurves{runCount} = result.yIntercepts;
      waitbar(runCount/totalRuns, h3, ...
        sprintf('Grow/Bright: %d / %d (%.1f%%)', runCount, totalRuns, 100*runCount/totalRuns));
    end
  end
end
close(h3);

% save all curves to disk
save('allExperiments_growth‐vs‐brightnessCurves.mat', 'allCurves');

%% 4) ROI‐Spacing Sweep
roiConfigs = {
  4, 2.0, 2.5;
  6, 2.5, 3.0;
  8, 2.0, 3.5
};
nConfigs  = size(roiConfigs,1);
totalRuns = nConfigs * runsPerCfg;
h4 = waitbar(0, 'ROI‐spacing sweep: 0%');
baseParams.testName         = 'ROI‐Spacing';



runCount = 0;
for idx = 1:nConfigs
  Nroi = roiConfigs{idx,1};
  a    = roiConfigs{idx,2};
  b    = roiConfigs{idx,3};
  for iter = 1:runsPerCfg
    runCount = runCount + 1;
    params = baseParams;
    params.radiiFactors = linspace(a,b,Nroi);
    params.seed         = uint32(Nroi*1e3 + round(a*10)*10 + iter);
    runDynamicTIS(params, false);
    waitbar(runCount/totalRuns, h4, ...
      sprintf('ROI spacing: %d / %d (%.1f%%)', runCount, totalRuns, 100*runCount/totalRuns));
  end
end
close(h4);

% save all curves to disk
save('allExperiments_roi‐spacingCurves.mat', 'allCurves');

