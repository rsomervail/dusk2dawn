% 
% 
%       plot validation metrics from up to 3 parameter dimensions at a time
%           inputs:
%               - cfg.dims2plot = indexes of which parameters will be plotted and how (1st is x-axis, 2nd is subplot rows, 3rd is subplot columns)
%           - ! make sure all figures match the single subject plot, 
%           - ! e.g. add the combined computation plot to single sub function
%           - ! check that the subplot indexing works properly with 3dims (i.e. that it doesn't transpose or something weird)
%           - ! write documentation
% 
%               - Richard Somervail, 22/11/2021 
%                   www.iannettilab.net
function d2d_groupLevel_plotValidation(cfg)

%% defaults
% general
if ~isfield(cfg,'loadPath'),        cfg.loadPath  = [];     end % set later
if ~isfield(cfg,'saveFig'),         cfg.saveFig   = true;   end
if ~isfield(cfg,'saveAs'),          cfg.saveAs    = '';     end
if ~isfield(cfg,'saveFold'),        cfg.saveFold  = [];     end % set later
if ~isfield(cfg,'maximiseFigs'),    cfg.maximiseFigs  = true;     end 

% what to plot
if ~isfield(cfg,'plot_asr'),    cfg.plot_asr  = true; end
if ~isfield(cfg,'plot_fft'),    cfg.plot_fft  = true; end
if ~isfield(cfg,'plot_sw'),     cfg.plot_sw   = true; end
if ~isfield(cfg,'plot_pica'),   cfg.plot_pica = true; end
if ~isfield(cfg,'plot_ica'),    cfg.plot_ica  = true; end
% which dimensions to plot and which ranges?
if ~isfield(cfg,'dims2plot'), cfg.dims2plot = []; end % set later
if ~isfield(cfg,'plotRange'), cfg.plotRange = []; end % set later

%% constants 
lms = 10;
MarkerSize = 4;
LineWidth  = 1.5;
AreaAlpha  = 0.5; 
SingleSubOpacity = 0.4;

cols_classProb4 = [ 0,0,1; 1,0.5,0; 0,0,0; 1,0,0  ]; % brain=blue, bioArtifact=orange, recArtifact=black, unknown=red

%% prepare stuff
% check for validation structure containing group-level validation metrics
if isempty(cfg.loadPath) && ~isfield(cfg,'validG')
    [file, path] = uigetfile(cd,'choose file containing group-level validation metrics');
    load([path '\' file], 'validG');
    g = validG;
elseif ~isempty(cfg.loadPath) && ~isfield(cfg,'validG')
    load([cfg.loadPath], 'validG');
    g = validG;
elseif isfield(cfg,'validG')
    g = cfg.validG;
end
npars  = length(g.pars.labels);
nfiles = length(g.datasets); 

% set default dimensions to plot
if isempty(cfg.dims2plot), cfg.dims2plot = 1:npars; end
if isempty(cfg.plotRange), cfg.plotRange = repmat({':'},1,npars); end

% get x-axis cleaning parameter
dimx = cfg.dims2plot(1);
parx_lab = g.pars.labels{dimx};
parx_lab = strrep(parx_lab,'_','-'); % replace underscores for plotting
parx_val = cellfun(@num2str, g.pars.values{dimx}, 'UniformOutput', false); % formatted as cell string array so that it can be applied to the labels regardless of numerical dimensions
parx_val = strrep(parx_val, 'NaN','raw'); % replace NaNs with more helpful label for original data
parx_val = regexprep(parx_val, '\s*',' '); % trim excess spaces#

% get 2nd and 3rd cleaning parameter dimensions (for subplots)
if length( cfg.dims2plot ) >= 2
    dimc = cfg.dims2plot(2); 
    parc_lab = g.pars.labels{dimc};
    parc_lab = strrep(parc_lab,'_','-'); % replace underscores for plotting
    parc_val = cellfun(@num2str, g.pars.values{dimc}, 'UniformOutput', false);
%     parc_val(1) = []; % remove original data value (redundant with x-axis plot)
    parc_val = regexprep(parc_val, '\s*',' '); % trim excess spaces
else
    dimc     = nan;
    parc_lab = nan;
    parc_val = nan;
end
if length( cfg.dims2plot ) >= 3
    dimr = cfg.dims2plot(3);
    parr_lab = g.pars.labels{dimr};
    parr_lab = strrep(parr_lab,'_','-'); % replace underscores for plotting
    parr_val = cellfun(@num2str, g.pars.values{dimr}, 'UniformOutput', false);
%     parr_val(1) = []; % remove original data value (redundant with x-axis plot)
    parr_val = regexprep(parr_val, '\s*',' '); % trim excess spaces
else
    dimr = nan;
    parr_lab = nan;
    parr_val = nan;
end

% get number of values spanning the 2 subplot dimensions
nparc = length(parc_val);
nparr = length(parr_val);
nparx = length(parx_val);

%% Initialise all figures

% cleanData
if cfg.plot_asr
    fig_ASRresults       = figure('name',[ g.name  ': how much did variance ASR remove?'], 'numbertitle','off');
    fig_computationTime  = figure('name',[ g.name  ': computation times'], 'numbertitle','off');
    fig_calibData        = figure('name',[ g.name  ': how much data did ASR have to calibrate with?'], 'numbertitle','off');
end
% fft
if isfield(g,'fft') && cfg.plot_fft
    fig_fft_powerByBand_abs  = figure('name',[ g.name  ': FFT power by band (absolute)'], 'numbertitle','off');
    fig_fft_powerByBand_prop = figure('name',[ g.name  ': FFT power by band (%)'], 'numbertitle','off');
end
% sw
if isfield(g,'sw')  && cfg.plot_sw
    fig_sw_nwaves         =  figure('name',[ g.name  ': SW number of slow-waves'], 'numbertitle','off');
    fig_sw_nwaves_verified_prop  =  figure('name',[ g.name  ': SW % of slow-waves verified'], 'numbertitle','off');
    fig_sw_mean_amp_abs   =  figure('name',[ g.name  ': SW slow-wave amplitude (absolute)'], 'numbertitle','off');
    fig_sw_mean_amp_prop  =  figure('name',[ g.name  ': SW slow-wave amplitude (%)'], 'numbertitle','off');
    fig_sw_tval_median    =  figure('name',[ g.name  ': SW slow-wave consistency (median t value)'], 'numbertitle','off');
end
if isfield(g,'pica') && cfg.plot_pica
    fig_pica_numIC                      =  figure('name',[ g.name  ': PICA number of components'], 'numbertitle','off');
%     fig_pica_tElapsed_abs               =  figure('name',[ g.name  ': PICA computation time (absolute)'], 'numbertitle','off');
    fig_pica_tElapsed_prop              =  figure('name',[ g.name  ': PICA computation time (percentage)'], 'numbertitle','off');
    fig_pica_identiDiff_meandiff        =  figure('name',[ g.name  ': PICA identifiability score'], 'numbertitle','off');
    fig_pica_identiDiff_propNumKnown    =  figure('name',[ g.name  ': PICA % ICs which were identifiable'], 'numbertitle','off');
    fig_pica_identiDiff_propVarKnown    =  figure('name',[ g.name  ': PICA % variance explained by identifiable ICs'], 'numbertitle','off');
    fig_pica_classProb4_catNumProp      =  figure('name',[ g.name  ': PICA % ICs of each category'], 'numbertitle','off');
    fig_pica_classProb4_catVarProp      =  figure('name',[ g.name  ': PICA % variance explained by each category'], 'numbertitle','off');
end
if isfield(g,'ica') && cfg.plot_ica
%     fig_ica_tElapsed_abs              =  figure('name',[ g.name   ': ICA computation time (absolute)'], 'numbertitle','off');
    fig_ica_tElapsed_prop             =  figure('name',[ g.name   ': ICA computation time (percentage)'], 'numbertitle','off');
    fig_ica_identiDiff_meandiff        =  figure('name',[ g.name  ': ICA identifiability score'], 'numbertitle','off');
    fig_ica_identiDiff_propNumKnown    =  figure('name',[ g.name  ': ICA % ICs which were identifiable'], 'numbertitle','off');
    fig_ica_identiDiff_propVarKnown    =  figure('name',[ g.name  ': ICA % variance explained by identifiable ICs'], 'numbertitle','off');
    fig_ica_classProb4_catNumProp      =  figure('name',[ g.name  ': ICA % ICs of each category'], 'numbertitle','off');
    fig_ica_classProb4_catVarProp      =  figure('name',[ g.name  ': ICA % variance explained by each category'], 'numbertitle','off');
end

%% Loop through cleaning parameters which are plotted on seperate figures
% get range to loop through for the row and column dimensions (x axis is handled below)
if isnan(dimr) 
    dimr_range = 1; % if not even plotting this dimension, just loop through 1
elseif strcmp(cfg.plotRange{dimr},':')
    dimr_range = 1:nparr; % use entire range of values
else
    eval([ 'dimr_range = ' cfg.plotRange{dimc} ';' ]) % use subset of values  
end
if isnan(dimc) 
    dimc_range = 1; % if not even plotting this dimension, just loop through 1
elseif strcmp(cfg.plotRange{dimc},':')
    dimc_range = 1:nparc; % use entire range of values
else
    eval([ 'dimc_range = ' cfg.plotRange{dimc} ';' ]) % use subset of values  
end

% get x axis parameter value range to plot
if strcmp(cfg.plotRange{dimx},':')
    dimx_range = 1:nparx; % use entire range of values
else
    eval([ 'dimx_range = ' cfg.plotRange{dimx} ';' ]) % use subset of values  
end

% change nparr, nparc, nparx to refer to the length of the parameter values actually being plotted
nparr = length(dimr_range);
nparc = length(dimc_range);
% nparx = length(dimx_range); % ? actually not used but written just in case

for r = 1:nparr
    for c = 1:nparc
        
%% construct a plot-index string for accessing validation matrices specific for this iteration of the row/column loop
indstr = cfg.plotRange; % get plot range for x axis (row and column ranges will be overwritten by the values for this iteration)
if ~isnan(dimc), indstr{dimc} = num2str(dimc_range(c)+1); end % after these two steps, should be left with one index corresponding to the x axis range
if ~isnan(dimr), indstr{dimr} = num2str(dimr_range(r)+1); end % +1 because first element of each dim is always original data
% if isnan(dimr) && isnan(dimc), indstr{2} = '2'; end % !!! causes errors??? ! if only using x dimension, need this extra index to ignore the "uncleaned only" column
indstr = [':,' strjoin(indstr,',')]; % make into string, first colon refers to dimension of merged files

%% construct title for each subplot
switch length( cfg.dims2plot )
    case 1
        subplot_title = [];
    case 2
        subplot_title = [parc_lab ': ' parc_val{dimc_range(c)+1}  ];
    case 3
        subplot_title = [parc_lab ': ' parc_val{dimc_range(c)+1} '\t' parr_lab ': ' parr_val{dimr_range(r)+1}   ];
end

%% PLOT - ASR cleaning metrics
if isfield(g,'asr') && cfg.plot_asr
    o = g.asr; % get subfield

    % percentage of timepoints altered by ASR
    set(0, 'CurrentFigure', fig_ASRresults) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    yyaxis left
    eval([ 'temp2plot = o.g_propCleaned(' indstr ');' ]) 
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
        axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel '% timepoints altered'
    ylim([0 100])

    % percentage of variance removed by ASR
    yyaxis right
    eval([ 'temp2plot = o.g_propVarRemoved(' indstr ');' ])
%     eval([ 'temp2plot = o.g_propGFPRemoved(' indstr ');' ])
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o', 'cmap',[1,0.4,0]  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    ylabel '% variance removed'
    ylim([0 100])
    title(subplot_title);
    
    % computation times - ASR part
    set(0, 'CurrentFigure', fig_computationTime) % select this figure
    eval([ 'temp2plot = o.g_tElapsed(' indstr ');' ])
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    temp2plot = double(temp2plot);
    temp2plot = temp2plot ./ 60; % change units to minutes
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel 'computation time (mins)'
    title(subplot_title);
    
    
    % how sufficient was calibration data?
    % minimum calibration times
    set(0, 'CurrentFigure', fig_calibData) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    yyaxis left
    eval([ 'temp2plot = o.g_calibDataLen_min(' indstr ');' ]) 
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
        axis square
    for sb = 1:size(temp2plot,1), plot( 1:length(dimx_range), temp2plot(sb,:), '-', 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel 'minimum calibration time (s)'
    ylim([ min(o.g_calibDataLen_min,[],'all')  max(o.g_calibDataLen_min,[],'all')])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
    title(subplot_title);
    plot(xlim,[60,60],'r--'); % 60s threshold dashed red line
    % percentage of chunks with adequate calibration data
    yyaxis right
    eval([ 'temp2plot = o.g_calibDataLen_prop60s(' indstr ');' ])
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o', 'cmap',[1,0.4,0]  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    ylabel '% of cleaning chunks with sufficient calibration data (>=60s)'
    ylim([0 100])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
    title(subplot_title);

    drawnow
end
%% PLOT - frequency decomposition metrics
if isfield(g,'fft') && cfg.plot_fft
    o = g.fft; % get subfield

    % plot frequency power by band (absolute)
    set(0, 'CurrentFigure', fig_fft_powerByBand_abs) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = squeeze( o.g_binAmp_avgChan(' indstr ',:) );' ])
    temp2plot = double(temp2plot);
    cols = distinguishable_colors(size(o.binFreqs,1));
    clear lineh areah
    for k = 1:size(temp2plot,3)
        [lineh(k), areah(k)] = boundedline( 1:length(dimx_range), mean(temp2plot(:,:,k)), std(temp2plot(:,:,k)) ,'-o', 'cmap', cols(k,:)  ); hold on
            lineh(k).LineWidth  = LineWidth;
            lineh(k).MarkerSize = MarkerSize;
            lineh(k).MarkerFaceColor = lineh(k).Color;
            areah(k).FaceAlpha  = AreaAlpha;
    end
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel 'power'
    if c == 1 && r == 1, legend(lineh, o.binFreqsLabels, 'Location','northeast'); end
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 
    
    
    % plot frequency power by band (percentage)
    set(0, 'CurrentFigure', fig_fft_powerByBand_prop) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    title(subplot_title); hold on
    eval([ 'temp2plot = squeeze( o.g_binAmp_avgChan_prop(' indstr ',:) );' ])
    temp2plot = double(temp2plot);
    cols = distinguishable_colors(size(o.binFreqs,1));
    clear lineh areah
    for k = 1:size(temp2plot,3)
        [lineh(k), areah(k)] = boundedline( 1:length(dimx_range), mean(temp2plot(:,:,k)), std(temp2plot(:,:,k)) ,'-o', 'cmap', cols(k,:)  ); hold on
            lineh(k).LineWidth  = LineWidth;
            lineh(k).MarkerSize = MarkerSize;
            lineh(k).MarkerFaceColor = lineh(k).Color;
            areah(k).FaceAlpha  = AreaAlpha;
    end
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel 'power (%)'
    ylim([0 120])
    if c == 1 && r == 1, legend(lineh, o.binFreqsLabels, 'Location','southwest'); end
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 

    drawnow
end
%% PLOT - slow-wave metrics
if isfield(g,'sw') && cfg.plot_sw
    o = g.sw; % get subfield

    % plot n slow-waves
    set(0, 'CurrentFigure', fig_sw_nwaves) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = o.g_nwaves(' indstr ');' ])
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    for sb = 1:size(temp2plot,1), plot( 1:length(dimx_range), temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel 'number of slow-waves'
    ylim([ min(o.g_nwaves,[],'all')  max(o.g_nwaves,[],'all')])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms])
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 
    
    % plot n slow-waves verified (percentage)
    set(0, 'CurrentFigure', fig_sw_nwaves_verified_prop) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = o.g_nwaves_verified_prop(' indstr ');' ])
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    for sb = 1:size(temp2plot,1), plot( 1:length(dimx_range), temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel '% of slow-waves verified'
    ylim([ 0, 100])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 
    
    % plot slow-wave amplitude (absolute)
    set(0, 'CurrentFigure', fig_sw_mean_amp_abs) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = o.g_amp_mean(' indstr ');' ])
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    for sb = 1:size(temp2plot,1), plot( 1:length(dimx_range), temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel 'slow-wave amplitude (µV)'
    ylim([ min(o.g_amp_mean,[],'all')  max(o.g_amp_mean,[],'all')])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 
    
    % plot slow-wave amplitude (percentage)
    set(0, 'CurrentFigure', fig_sw_mean_amp_prop) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = o.g_amp_mean_prop(' indstr ');' ])
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    for sb = 1:size(temp2plot,1), plot( 1:length(dimx_range), temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel 'slow-wave amplitude (%)'
    ylim([ 0, 100])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 
    
    % plot slow-wave t-value
    set(0, 'CurrentFigure', fig_sw_tval_median) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = o.g_tval_median(' indstr ');' ])
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    for sb = 1:size(temp2plot,1), plot( 1:length(dimx_range), temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel 'slow-wave consistency (median t value)'
    ylim([ min(o.g_tval_median,[],'all')  max(o.g_tval_median,[],'all')])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms])
    set(gca, 'YDir', 'reverse')
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 

    drawnow
end
%% PLOT - pica metrics
if isfield(g,'pica') && cfg.plot_pica
    o = g.pica; % get subfield

    % plot number of PICA components
    set(0, 'CurrentFigure', fig_pica_numIC) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = o.g_numIC(' indstr ');' ])
    temp2plot = double(temp2plot);
    temp2plot = temp2plot ./ 60; % change units to minutes
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    for sb = 1:size(temp2plot,1), plot( 1:length(dimx_range), temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel 'number of ICs'
    ylim([ min(o.g_numIC,[],'all')  max(o.g_numIC,[],'all')])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 
    
    % plot computation time (absolute)
    set(0, 'CurrentFigure', fig_computationTime) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = o.g_tElapsed(' indstr ');' ])
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-og'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel 'computation time (mins)'
    title(subplot_title);
    
    % plot computation time (percentage)
    set(0, 'CurrentFigure', fig_pica_tElapsed_prop) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = o.g_tElapsed_prop(' indstr ');' ])
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    for sb = 1:size(temp2plot,1), plot( 1:length(dimx_range), temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel 'computation time (%)'
    ylim([0 100])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 
    
    % plot identifiability (meanDiff)
    set(0, 'CurrentFigure', fig_pica_identiDiff_meandiff) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = o.g_identiDiff_meanDiff(' indstr ');' ])
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    for sb = 1:size(temp2plot,1), plot( 1:length(dimx_range), temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel 'mean identifiability score'
    ylim([ -1  1])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 
    
    % plot identifiability (propNumKnown)
    set(0, 'CurrentFigure', fig_pica_identiDiff_propNumKnown) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = o.g_identiDiff_propNumKnown(' indstr ');' ])
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    for sb = 1:size(temp2plot,1), plot( 1:length(dimx_range), temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel '% identifiable components'
    ylim([ 0 100 ])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 
    
    % plot identifiability (propVarKnown)
    set(0, 'CurrentFigure', fig_pica_identiDiff_propVarKnown) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = o.g_identiDiff_propVarKnown(' indstr ');' ])
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    for sb = 1:size(temp2plot,1), plot( 1:length(dimx_range), temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel '% variance explained by identifiable components'
    ylim([ 0 100 ])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 
    
    % plot class probabilities (% of components)
    set(0, 'CurrentFigure', fig_pica_classProb4_catNumProp) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = squeeze( o.g_classProb4_catNumProp(' indstr ',:) );' ])
    temp2plot = double(temp2plot);
    clear lineh areah
    for k = 1:size(temp2plot,3)
        [lineh(k), areah(k)] = boundedline( 1:length(dimx_range), mean(temp2plot(:,:,k)), std(temp2plot(:,:,k)) ,'-o', 'cmap', cols_classProb4(k,:)  ); hold on
            lineh(k).LineWidth  = LineWidth;
            lineh(k).MarkerSize = MarkerSize;
            areah(k).FaceAlpha  = AreaAlpha;
    end
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel '% of components'
    ylim([0 100])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
    if c == 1 && r == 1, legend(lineh, o.classProb4_classLabs, 'Location','northeast'); end
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 
    
    % plot class probabilities (% of variance explained)
    set(0, 'CurrentFigure', fig_pica_classProb4_catVarProp) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = squeeze( o.g_classProb4_catVarProp(' indstr ',:) );' ])
    temp2plot = double(temp2plot);
    clear lineh areah
    for k = 1:size(temp2plot,3)
        [lineh(k), areah(k)] = boundedline( 1:length(dimx_range), mean(temp2plot(:,:,k)), std(temp2plot(:,:,k)) ,'-o', 'cmap', cols_classProb4(k,:)  ); hold on
            lineh(k).LineWidth  = LineWidth;
            lineh(k).MarkerSize = MarkerSize;
            areah(k).FaceAlpha  = AreaAlpha;
    end
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel '% of variance explained'
    ylim([0 100])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
    if c == 1 && r == 1, legend(lineh, o.classProb4_classLabs, 'Location','northeast'); end
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 

    drawnow
end
%% PLOT - ica metrics
if isfield(g,'ica') && cfg.plot_ica
    o = g.ica; % get subfield

    % plot computation time (absolute)
    set(0, 'CurrentFigure', fig_computationTime) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = o.g_tElapsed(' indstr ');' ])
    temp2plot = double(temp2plot);
    temp2plot = temp2plot ./ 60; % change units to minutes
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-or'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel 'computation time (mins)'
    title(subplot_title);
    
    % plot computation time (percentage)
    set(0, 'CurrentFigure', fig_ica_tElapsed_prop) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = o.g_tElapsed_prop(' indstr ');' ])
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    for sb = 1:size(temp2plot,1), plot( 1:length(dimx_range), temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel 'computation time (%)'
    ylim([0 100])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 
    
    % plot identifiability (meanDiff)
    set(0, 'CurrentFigure', fig_ica_identiDiff_meandiff) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = o.g_identiDiff_meanDiff(' indstr ');' ])
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    for sb = 1:size(temp2plot,1), plot( 1:length(dimx_range), temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel 'mean identifiability score'
    ylim([ -1  1])
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 
    
    % plot identifiability (propNumKnown)
    set(0, 'CurrentFigure', fig_ica_identiDiff_propNumKnown) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = o.g_identiDiff_propNumKnown(' indstr ');' ])
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    for sb = 1:size(temp2plot,1), plot( 1:length(dimx_range), temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel '% identifiable components'
    ylim([ 0 100 ])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 
    
    % plot identifiability (propVarKnown)
    set(0, 'CurrentFigure', fig_ica_identiDiff_propVarKnown) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = o.g_identiDiff_propVarKnown(' indstr ');' ])
    temp2plot = double(temp2plot);
    [lineh, areah] = boundedline( 1:length(dimx_range), mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        lineh.LineWidth  = LineWidth;
        lineh.MarkerSize = MarkerSize;
        lineh.MarkerFaceColor = lineh.Color;
        areah.FaceAlpha  = AreaAlpha;
    for sb = 1:size(temp2plot,1), plot( 1:length(dimx_range), temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel '% variance explained by identifiable components'
    ylim([ 0 100 ])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 
    
    % plot class probabilities (% of components)
    set(0, 'CurrentFigure', fig_ica_classProb4_catNumProp) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = squeeze( o.g_classProb4_catNumProp(' indstr ',:) );' ])
    temp2plot = double(temp2plot);
    clear lineh areah
    for k = 1:size(temp2plot,3)
        [lineh(k), areah(k)] = boundedline( 1:length(dimx_range), mean(temp2plot(:,:,k)), std(temp2plot(:,:,k)) ,'-o', 'cmap', cols_classProb4(k,:)  ); hold on
            lineh(k).LineWidth  = LineWidth;
            lineh(k).MarkerSize = MarkerSize;
            lineh(k).MarkerFaceColor = lineh(k).Color;
            areah(k).FaceAlpha  = AreaAlpha;
    end
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel '% of components'
    ylim([0 100])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
    if c == 1 && r == 1, legend(lineh,o.classProb4_classLabs, 'Location','northeast'); end
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 
    
    % plot class probabilities (% of variance explained)
    set(0, 'CurrentFigure', fig_ica_classProb4_catVarProp) % select this figure
    subplot(nparr,nparc, getSubplotInd(r,c, nparc));
    eval([ 'temp2plot = squeeze( o.g_classProb4_catVarProp(' indstr ',:) );' ])
    temp2plot = double(temp2plot);
    clear lineh areah
    for k = 1:size(temp2plot,3)
        [lineh(k), areah(k)] = boundedline( 1:length(dimx_range), mean(temp2plot(:,:,k)), std(temp2plot(:,:,k)) ,'-o', 'cmap', cols_classProb4(k,:)  ); hold on
            lineh(k).LineWidth  = LineWidth;
            lineh(k).MarkerSize = MarkerSize;
            lineh(k).MarkerFaceColor = lineh(k).Color;
            areah(k).FaceAlpha  = AreaAlpha;
    end
    axis square
    xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylabel '% of variance explained'
    ylim([0 100])
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
    if c == 1 && r == 1, legend(lineh,o.classProb4_classLabs, 'Location','northeast'); end
    title(subplot_title);
    reorderPlots(gca); % bring all lines to foreground and send patches to background 

    drawnow
end

%% end loop through parameters plotted on seperate figures
    end % end row loop
end % end col loop

%% set legend for computation time figure
if any( [ cfg.plot_asr , cfg.plot_ica, cfg.plot_pica] )
    set(0, 'CurrentFigure', fig_computationTime) % select this figure
    compLegend = { 'ASR', 'PICA', 'ICA' };
    compLegend = compLegend( [ cfg.plot_asr , cfg.plot_pica, cfg.plot_ica] & [ isfield(g,'cleanData'), isfield(g,'pica'), isfield(g,'ica')]   );
    ax = fig_computationTime.Children(1);
    for p = 1:length(ax.Children), legPlots{p} = ax.Children(p).Type; end
    legPlots = ax.Children(  strcmp( legPlots , 'line' ) );
    legend(  legPlots(end:-1:1)  , compLegend , 'Location','northeast');
    axs = findall( fig_computationTime, 'type', 'axes'  );
    ylims = [min([axs.YLim]) max([axs.YLim])];
    for k = 1:length(axs)
       axs(k).YLim = ylims; 
    end
end

%% Maximise all figures
if cfg.maximiseFigs
    figs = findall(0, 'type', 'Figure');
    for f = 1:length(figs)
        fig = figs(f);
        fig.Units = 'normalized';
        fig.Position = [0 0.0370 1 0.8917];
    end
end

%% Save figures
if cfg.saveFig
    fprintf('saving figures ...\n')
    
    % get default savefolder if no folder is specified
    if isempty(cfg.saveFold)
        if exist('EEG','var')
            temp = EEG.etc.cleanSleep.valid.loadPath; % use location of uncleaned data if no path specified
            cfg.saveFold = temp(1:end-(length(EEG.filename)));
        else
            cfg.saveFold = cd; % set to cd if no path can be found from uncleaned data
        end
    end
    
    % loop through figure handles and save
    figs = who('fig_*'); % get all figure handle variable names
    temp = strsplit(cfg.loadPath,{'\','/'}); temp = temp{end}(1:end-4);
    for f = 1:length(figs)
        eval([ 'savefig(' figs{f}  ', [''fig_' temp '_'  cfg.saveAs '_'  figs{f}(5:end) '''])' ]) % save fig
    end
    
    fprintf('... saved figures\n')
end

%% subfunctions
function reorderPlots(ax)
    

    % assuming patch is always last
    lines = []; theRest = [];
    for j = 1:length(ax.Children)
        if strcmp( ax.Children(j).Type, 'line' )
            if ~strcmp( ax.Children(j).Marker , 'none' ) % average line
                groupline = j;
            else % single subject line
                lines = [lines, j];
            end
        else % everything else (patches for error bars) goes to the background
            theRest = [theRest, j];
        end
    end
    % change plot order
    ax.Children = ax.Children( [ groupline, lines, theRest ]  );

     
end


function ind = getSubplotInd(r,c,numc)
    ind = (r-1)*numc + c;
end


end % END MAIN FUNCTION

