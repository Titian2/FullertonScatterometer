function results = runDynamicTIS(params, doPlot)
% runDynamicTIS  Run the dynamic‐ROI toy model with arbitrary parameters.
%   results = runDynamicTIS(params, doPlot)
%
% Writes two TSVs:
%   <testName>_seed<seed>_summary.tsv
%   <testName>_seed<seed>_initialCircles.tsv

  arguments
    params struct
    doPlot (1,1) logical = true
  end

  %% 0) Seeding for reproducibility
  p = params;
  if isfield(p,'seed')
    rng(p.seed);
  else
    s = rng('shuffle');
    p.seed = s.Seed;
  end

  %% 1) Unpack & defaults
  if ~isfield(p,'embedSize'),   p.embedSize = 10; end
  Ncircles  = p.numRandomCircles;
  ITERS     = p.numIterations;
  SZ        = p.imgSize;
  IR        = p.initialRadius;
  RF        = p.radiiFactors;

  %% 2) Sample circle placement & amplitudes
  [Xbg,Ybg] = meshgrid(1:SZ,1:SZ);

  switch lower(p.placementBias.type)
    case 'uniform'
      cx = randi(SZ,[1,Ncircles]);
      cy = randi(SZ,[1,Ncircles]);
      biasParam = NaN;
    case 'gaussian'
      sigma = p.placementBias.sigma;
      cx = round(SZ/2 + sigma*randn(1,Ncircles));
      cy = round(SZ/2 + sigma*randn(1,Ncircles));
      cx = min(max(cx,1),SZ);
      cy = min(max(cy,1),SZ);
      biasParam = sigma;
    otherwise
      error('Unknown placementBias.type');
  end

  switch lower(p.initAmpDist.type)
    case 'uniform'
      a0 = p.initAmpDist.min + ...
           (p.initAmpDist.max-p.initAmpDist.min).*rand(1,Ncircles);
      ampParam1 = p.initAmpDist.min;
      ampParam2 = p.initAmpDist.max;
    case 'normal'
      a0 = p.initAmpDist.mu + p.initAmpDist.sigma.*randn(1,Ncircles);
      a0(a0<0)=0;
      ampParam1 = p.initAmpDist.mu;
      ampParam2 = p.initAmpDist.sigma;
    otherwise
      error('Unknown initAmpDist.type');
  end

  %% 3) Build initial image
  background = ones(SZ);
  for k=1:Ncircles
    mask = (Xbg-cx(k)).^2 + (Ybg-cy(k)).^2 <= IR^2;
    background(mask) = background(mask) + a0(k);
  end
  % central embed
  E = p.embedSize;
  patch = ones(E)*10 + 1;
  y0 = round(SZ/2 - E/2)+1;
  x0 = round(SZ/2 - E/2)+1;
  background(y0:y0+E-1, x0:x0+E-1) = patch;

  img = background;

  %% 4) Define ROIs
  width_s = p.embedSize * 1.1;
  radii   = width_s * RF;
  Nroi    = numel(radii);
  [X,Y]   = meshgrid(1:SZ,1:SZ);

  %% 5) Prep storage & optional plotting
  yIntercepts = nan(1,ITERS);

  if doPlot
    figure('Name',p.testName,'Position',[200 200 900 400]);
    ax1 = subplot(1,2,1); 
      hImg = imagesc(img,'Parent',ax1); colormap(ax1,'gray');
      axis(ax1,'image','off'); hold(ax1,'on');
      th = linspace(0,2*pi,360);
      cols = {'b--','r--','g--','y--','m--','c--'};
      for i=1:Nroi
        xc = SZ/2 + radii(i)*cos(th);
        yc = SZ/2 + radii(i)*sin(th);
        plot(ax1,xc,yc,cols{mod(i-1,6)+1},'LineWidth',2);
      end
      title(ax1,'ROIs','FontSize',14);

    ax2 = subplot(1,2,2);
      hPlot = plot(ax2,nan,nan,'b-','LineWidth',1.5);
      xlabel(ax2,'Iteration'); ylabel(ax2,'Y-Intercept');
      title(ax2,'Intercept vs Iter','FontSize',14);
      xlim(ax2,[1 ITERS]); grid(ax2,'on'); hold(ax2,'on');
  end

  %% 6) Main evolution loop
  for it = 1:ITERS
    % grow & brighten
    maskAll = false(SZ);
    for k=1:Ncircles
      rNow = IR + p.growthRate*it;
      maskAll = maskAll | ((Xbg-cx(k)).^2 + (Ybg-cy(k)).^2 <= rNow^2);
    end
    img(maskAll) = img(maskAll) + p.brightnessStep;

    % sums & areas
    sums = zeros(1,Nroi); areas = sums;
    for i=1:Nroi
      mROI    = (X-SZ/2).^2 + (Y-SZ/2).^2 <= radii(i).^2;
      sums(i) = sum(img(mROI));
      areas(i)= sum(mROI(:));
    end

    % fit and record
    c = polyfit(areas,sums,1);
    yIntercepts(it) = c(2);

    if doPlot
      set(hImg, 'CData', img);
      set(hPlot,'XData',1:it,'YData',yIntercepts(1:it));
      drawnow;
    end
  end

  %% 7) Diagnostics: circle‐in/out/boundary counts
  d.dist = sqrt((cx-SZ/2).^2 + (cy-SZ/2).^2);
  tol    = 0.5;
  d.insideSmall  = sum(d.dist < radii(1));
  d.insideLarge  = sum(d.dist < radii(end));
  Dmat           = abs(bsxfun(@minus,d.dist(:),radii(:)'));
  d.onBoundary   = sum(any(Dmat<=tol,2));
  d.outsideAll   = sum(d.dist > radii(end));

  %% 8) Package results
  results.params        = p;
  results.yIntercepts   = yIntercepts;
  results.metrics.seed              = p.seed;
  results.metrics.initialIntercept  = yIntercepts(1);
  results.metrics.finalIntercept    = yIntercepts(end);
  results.metrics.numCircles        = Ncircles;
  results.metrics.initialRadius     = IR;
  results.metrics.endRadius         = IR + p.growthRate*ITERS;
  results.metrics.placementBiasType = p.placementBias.type;
  results.metrics.placementBiasSigma= biasParam;
  results.metrics.initAmpType       = p.initAmpDist.type;
  results.metrics.initAmpParam1     = ampParam1;
  results.metrics.initAmpParam2     = ampParam2;
  results.metrics.growthRate        = p.growthRate;
  results.metrics.brightnessStep    = p.brightnessStep;
  results.metrics.ROISpacing        = diff(radii);
  results.metrics.insideSmallROI    = d.insideSmall;
  results.metrics.insideLargeROI    = d.insideLarge;
  results.metrics.onROIBoundaries   = d.onBoundary;
  results.metrics.outsideROIs       = d.outsideAll;
  results.metrics.isNegative = ...
    (results.metrics.finalIntercept < results.metrics.initialIntercept);

 
  %% 9) Dump tab-delimited summary (append to one file, real tabs)
  summaryFile = sprintf('%s_summary.tsv', p.testName);

  % Build header with actual tab characters and a newline
  headerStr = sprintf([ ...
    'testName\tseed\tinitInt\tfinalInt\t', ...
    'nCircles\trStart\trEnd\tbiasType\tbiasSigma\t', ...
    'ampType\tampParam1\tampParam2\t', ...
    'growRate\tbrightStep\tROIspacing\t', ...
    'inSmall\tinLarge\tonBound\toutAll\tisNegative\n' ...
  ]);

  % Build the data line, also using \t between fields
  dataLine = sprintf([ ...
    '%s\t%d\t%.6g\t%.6g\t%d\t%.6g\t%.6g\t%s\t%.6g\t%s\t', ...
    '%.6g\t%.6g\t%.6g\t%.6g\t"%s"\t%d\t%d\t%d\t%d\t%d\n' ], ...
    p.testName, ...
    p.seed, ...
    results.metrics.initialIntercept, ...
    results.metrics.finalIntercept, ...
    Ncircles, ...
    IR, ...
    results.metrics.endRadius, ...
    p.placementBias.type, ...
    biasParam, ...
    p.initAmpDist.type, ...
    ampParam1, ...
    ampParam2, ...
    p.growthRate, ...
    p.brightnessStep, ...
    strjoin(string(results.metrics.ROISpacing),','), ...
    results.metrics.insideSmallROI, ...
    results.metrics.insideLargeROI, ...
    results.metrics.onROIBoundaries, ...
    results.metrics.outsideROIs, ...
    double(results.metrics.isNegative) ...
  );

  % If the file doesn't already exist, write the header first
  if ~isfile(summaryFile)
    fid = fopen(summaryFile,'w');
    fprintf(fid, '%s', headerStr);
  else
    fid = fopen(summaryFile,'a');
  end

  % Append the new data line
  fprintf(fid, '%s', dataLine);
  fclose(fid);

  %% 10) Dump tab-delimited initial‐circles
  initFile = sprintf('%s_seed%d_initialCircles.tsv',p.testName,p.seed);
  fid2 = fopen(initFile,'w');
  fprintf(fid2,'cx\tcy\tr0\tinitAmp\n');
  for k=1:Ncircles
    fprintf(fid2,'%d\t%d\t%.6g\t%.6g\n', cx(k), cy(k), IR, a0(k));
  end
  fclose(fid2);

end