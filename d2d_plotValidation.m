% 
% 
%       plot validation metrics from up to 3 parameter dimensions at a time
%           inputs:
%               - cfg.dims2plot = indexes of which parameters will be plotted and how (1st is x-axis, 2nd is subplot rows, 3rd is subplot columns)
% 
%
% 
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net
% History:
% 19/01/2023 ver 1.0.0 Created
% 
%%  
function d2d_plotValidation(EEG, cfg)

fprintf('d2d_plotValidation: plotting effects of ASR on various validation metrics for dataset: %s\n',EEG.etc.dusk2dawn.cfg.origFile)

%% defaults
% general
if ~isfield(cfg,'saveFig'),     cfg.saveFig   = false;   end
if ~isfield(cfg,'saveAs'),      cfg.saveAs    = '';     end
if ~isfield(cfg,'saveFold'),    cfg.saveFold  = [];     end % set later
% what to plot (? currently unused)
if ~isfield(cfg,'plot_asr'),    cfg.plot_asr  = true; end
if ~isfield(cfg,'plot_fft'),    cfg.plot_fft  = true; end
if ~isfield(cfg,'plot_sw'),     cfg.plot_sw   = true; end
if ~isfield(cfg,'plot_ica'),    cfg.plot_ica  = true; end
% which dimensions to plot and which ranges?
if ~isfield(cfg,'dims2plot'), cfg.dims2plot = []; end % set later
if ~isfield(cfg,'plotRange'), cfg.plotRange = []; end % set later

%% constants
lms = 10;
MarkerSize = 6;

%% get dataset info
origFile = strrep(EEG.etc.dusk2dawn.cfg.origFile,'.set','');
valid_merged = EEG.etc.dusk2dawn.valid_merged;

%% get parameters to plot on each axis
pars = EEG.etc.dusk2dawn.valid_merged(1).pars;
npars = length(pars.labels);

% if no varied parameters then create virtual parameter and set this to plot on X axis
if npars == 0
    pars.labels = {'ASR'};
    pars.values{1} = {'raw','cleaned'};
    cfg.dims2plot = 1;
    cfg.plotRange = {':'};
end

% default plot range
if isempty(cfg.plotRange), cfg.plotRange = repmat({':'},1,npars); end

% arrange dimensions to plot in more convenient way
if ~isfield(cfg, 'dims2plot')
    cfg.dims2plot = [ cfg.plotX cfg.plotRows cfg.plotCols  ];
    cfg.dims2plot = cfg.dims2plot( [cfg.dims2plot]~=0 ); % remove 0 dims (which are not plotted)
end

% get x-axis cleaning parameter
dimx = cfg.dims2plot(1);
parx_lab = pars.labels{dimx};
parx_lab = strrep(parx_lab,'_','-'); % replace underscores for plotting
parx_val = d2d_getParList(pars.values{dimx}); % extract varied parameter values
parx_val = cellfun(@num2str, parx_val, 'UniformOutput', false); % formatted as cell string array so that it can be applied to the labels regardless of numerical dimensions
parx_val = strrep(parx_val, 'NaN','raw'); % replace NaNs with more helpful label for raw data
parx_val = regexprep(parx_val, '\s*',' '); % trim excess spaces#

% get 2nd and 3rd cleaning parameter dimensions (for subplots)
if length( cfg.dims2plot ) >= 2
    dimc = cfg.dims2plot(2); 
    parc_lab = pars.labels{dimc};
    parc_lab = strrep(parc_lab,'_','-'); % replace underscores for plotting
    parc_val = d2d_getParList(pars.values{dimc}); % extract varied parameter values
    parc_val = cellfun(@num2str, parc_val, 'UniformOutput', false);
    parc_val(1) = []; % remove original data value (redundant with x-axis plot)
    parc_val = regexprep(parc_val, '\s*',' '); % trim excess spaces
else
    dimc     = nan;
    parc_lab = nan;
    parc_val = nan;
end
if length( cfg.dims2plot ) >= 3
    dimr = cfg.dims2plot(3);
    parr_lab = pars.labels{dimr};
    parr_lab = strrep(parr_lab,'_','-'); % replace underscores for plotting
    parr_val = d2d_getParList(pars.values{dimr}); % extract varied parameter values
    parr_val = cellfun(@num2str, parr_val, 'UniformOutput', false);
    parr_val(1) = []; % remove original data value (redundant with x-axis plot)
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

%% Loop through any selected sleep stages
if ~isfield(cfg,'stages')
    cfg.stages = 1;
    nstages = 1;
else
    stageCodes = [EEG.etc.dusk2dawn.valid_merged(:).stage];
    nstages = length(stageCodes);
end

for st = cfg.stages

    %% Initialise all figures

    v = valid_merged(st);
    if nstages > 1
        figprefix = [ origFile '_' stageCodes{st} ': ' ];
    else
        figprefix = [ origFile ': ' ];
    end

    % cleanData
    if cfg.plot_asr
        fig_cleanData(st)       = figure('name',[ figprefix 'how much variance did ASR remove?'], 'numbertitle','off');
        fig_calibData(st)       = figure('name',[ figprefix  'how much data did ASR have to calibrate with?'], 'numbertitle','off');
    end
    % fft
    if isfield(v,'fft') && cfg.plot_fft
        fig_fft_powerByBand_abs(st)  = figure('name',[ figprefix  'FFT power by band (absolute)'], 'numbertitle','off');
        fig_fft_powerByBand_prop(st) = figure('name',[ figprefix  'FFT power by band (%)'], 'numbertitle','off');
    end
    % sw
    if isfield(v,'sw')  && cfg.plot_sw
        fig_sw_nwaves(st)                   =  figure('name',[ figprefix  'SW number of slow-waves'], 'numbertitle','off');
        fig_sw_amp_mean_abs(st)     =  figure('name',[ figprefix  'SW slow-wave amplitude (absolute)'], 'numbertitle','off');
        fig_sw_amp_mean_prop(st)    =  figure('name',[ figprefix  'SW slow-wave amplitude (%)'], 'numbertitle','off');
        fig_sw_tval_median(st)      =  figure('name',[ figprefix  'SW slow-wave median t value'], 'numbertitle','off');
    end
    if isfield(v,'ica') && cfg.plot_ica
        fig_ica_tElapsed_abs(st)              =  figure('name',[ figprefix  'ICA computation time (absolute)'], 'numbertitle','off');
        fig_ica_tElapsed_prop(st)             =  figure('name',[ figprefix  'ICA computation time (percentage)'], 'numbertitle','off');
%         fig_ica_identiDiff_meandiff(st)        =  figure('name',[ figprefix  'ICA identifiability score'], 'numbertitle','off');
%         fig_ica_identiDiff_propNumKnown(st)    =  figure('name',[ figprefix  'ICA % ICs which were identifiable'], 'numbertitle','off');
%         fig_ica_identiDiff_propVarKnown(st)    =  figure('name',[ figprefix  'ICA % variance explained by identifiable ICs'], 'numbertitle','off');
        fig_ica_classProb4_catNumProp(st)      =  figure('name',[ figprefix  'ICA % ICs of each category'], 'numbertitle','off');
        fig_ica_classProb4_catVarProp(st)      =  figure('name',[ figprefix  'ICA % variance explained by each category'], 'numbertitle','off');
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
    indstr = strjoin(indstr,','); % make into string
    
    %% construct title for each subplot
    switch length( cfg.dims2plot )
        case 1
            subplot_title = [];
        case 2
            subplot_title = [parc_lab ': ' parc_val{dimc_range(c)}  ];
        case 3
            subplot_title = [parc_lab ': ' parc_val{dimc_range(c)} '\t' parr_lab ': ' parr_val{dimr_range(r)+1}   ];
    end
    
    %% PLOT - ASR cleaning metrics
    if isfield(v,'asr') && cfg.plot_asr
        o = v.asr; % get subfield
    
        % percentage of timepoints altered by ASR
        set(0, 'CurrentFigure', fig_cleanData(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        yyaxis left
        eval([ 'temp2plot = o.m_propCleaned(' indstr ');' ])
        plot( 1:length(dimx_range) , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; % ? perhaps plot solid line only if number of thresholds exceeds 3?
        axis square
        xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
        xlabel(parx_lab) % set x axis label to match parameter name
        xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
        ylabel '% timepoints altered'
        ylim([-10 110])
        title(subplot_title);
    
        % percentage of variance removed by ASR
        yyaxis right
        eval([ 'temp2plot = o.m_propVarRemoved(' indstr ');' ])
    %     eval([ 'temp2plot = o.m_propGFPRemoved(' indstr ');' ])
        plot( 1:length(dimx_range) , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; 
        ylabel '% variance removed'
        ylim([-10 110])
%         ylim([ min(o.m_propVarRemoved,[],'all')  max(o.m_propVarRemoved,[],'all')])
%         ylim([ min(o.m_propGFPRemoved,[],'all')  max(o.m_propGFPRemoved,[],'all')])
%         ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
        title(subplot_title);
        
        % how sufficient was calibration data?
        % minimum calibration times
        set(0, 'CurrentFigure', fig_calibData(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        yyaxis left
        eval([ 'temp2plot = o.m_calibDataLen_min(' indstr ');' ])
        plot( 1:length(dimx_range) , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; % ? perhaps plot solid line only if number of thresholds exceeds 3?
        axis square
        xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
        xlabel(parx_lab) % set x axis label to match parameter name
        xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
        ylabel 'minimum calibration time (s)'
        ylim([ min(o.m_calibDataLen_min,[],'all')  max(o.m_calibDataLen_min,[],'all')])
        ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
        title(subplot_title);
        plot(xlim,[60,60],'r--'); % 60s threshold dashed red line
        % percentage of data chunks with sufficient calibration
        yyaxis right
        eval([ 'temp2plot = o.m_calibDataLen_prop60s(' indstr ');' ])
        plot( 1:length(dimx_range) , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; 
        ylabel '% of cleaning chunks with sufficient calibration data (>=60s)'
        ylim([0 100])
        ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
        title(subplot_title);
    
        drawnow
    end
    %% PLOT - frequency decomposition metrics
    if isfield(v,'fft') && cfg.plot_fft
        o = v.fft; % get subfield
    
        % plot frequency power by band (absolute)
        set(0, 'CurrentFigure', fig_fft_powerByBand_abs(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        eval([ 'temp2plot = squeeze( o.m_binAmp_avgChan(' indstr ',:) );' ])
        for k = 1:size(temp2plot,2)
            plot( 1:length(dimx_range) , temp2plot(:,k), '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; % ? perhaps plot solid line only if number of thresholds exceeds 3?
        end
        axis square
        xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
        xlabel(parx_lab) % set x axis label to match parameter name
        xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
        ylabel 'power'
        if c == 1 && r == 1, legend(o.binFreqsLabels, 'Location','northeast'); end
        title(subplot_title);
        
        % plot frequency power by band (percentage)
        set(0, 'CurrentFigure', fig_fft_powerByBand_prop(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        title(subplot_title); hold on
        eval([ 'temp2plot = squeeze( o.m_binAmp_avgChan_prop(' indstr ',:) );' ])
        for k = 1:size(temp2plot,2)
            plot( 1:length(dimx_range) , temp2plot(:,k), '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; % ? perhaps plot solid line only if number of thresholds exceeds 3?
        end
        axis square
        xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
        xlabel(parx_lab) % set x axis label to match parameter name
        xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
        ylabel 'power (%)'
        ylim([ min(o.m_binAmp_avgChan_prop,[],'all')  max(o.m_binAmp_avgChan_prop,[],'all')])
        ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
        if c == 1 && r == 1, legend(o.binFreqsLabels, 'Location','northeast'); end
        title(subplot_title);
       
        drawnow
    end
    %% PLOT - slow-wave metrics
    if isfield(v,'sw') && cfg.plot_sw
        o = v.sw; % get subfield
    
        % plot n slow-waves
        set(0, 'CurrentFigure', fig_sw_nwaves(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        eval([ 'temp2plot = o.m_nwaves(' indstr ');' ])
        plot( 1:length(dimx_range) , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; % ? perhaps plot solid line only if number of thresholds exceeds 3?
        axis square
        xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
        xlabel(parx_lab) % set x axis label to match parameter name
        xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
        ylabel 'number of slow-waves'
        if (std(std(o.m_nwaves))~=0), ylim([ min(o.m_nwaves,[],'all')  max(o.m_nwaves,[],'all')]); ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]); end
        title(subplot_title);
        
        % plot slow-wave amplitude (absolute)
        set(0, 'CurrentFigure', fig_sw_amp_mean_abs(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        eval([ 'temp2plot = o.m_amp_mean(' indstr ');' ])
        plot( 1:length(dimx_range) , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; % ? perhaps plot solid line only if number of thresholds exceeds 3?
        axis square
        xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
        xlabel(parx_lab) % set x axis label to match parameter name
        xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
        ylabel 'slow-wave amplitude (µV)'
        ylim([ min(o.m_amp_mean,[],'all')  max(o.m_amp_mean,[],'all')])
        ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
        title(subplot_title);
        
        % plot slow-wave amplitude (percentage)
        set(0, 'CurrentFigure', fig_sw_amp_mean_prop(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        eval([ 'temp2plot = o.m_amp_mean_prop(' indstr ');' ])
        plot( 1:length(dimx_range) , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; % ? perhaps plot solid line only if number of thresholds exceeds 3?
        axis square
        xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
        xlabel(parx_lab) % set x axis label to match parameter name
        xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
        ylabel 'slow-wave amplitude (% of raw)'
        ylim([ 0, 100])
        ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
        title(subplot_title);
        
        % plot slow-wave t-value
        set(0, 'CurrentFigure', fig_sw_tval_median(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        eval([ 'temp2plot = o.m_tval_median(' indstr ');' ])
        plot( 1:length(dimx_range) , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; % ? perhaps plot solid line only if number of thresholds exceeds 3?
        axis square
        xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
        xlabel(parx_lab) % set x axis label to match parameter name
        xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
        ylabel 'slow-wave consistency (t value)'
        ylim([ min(o.m_tval_median,[],'all')  max(o.m_tval_median,[],'all')])
        ylim(ylim + [-range(ylim)/lms,range(ylim)/lms])
        title(subplot_title);
    
        drawnow
    end
    
    %% PLOT - ica metrics
    if isfield(v,'ica') && cfg.plot_ica
        o = v.ica; % get subfield
    
        % plot computation time (absolute)
        set(0, 'CurrentFigure', fig_ica_tElapsed_abs(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        eval([ 'temp2plot = o.m_tElapsed(' indstr ');' ])
        temp2plot = temp2plot / 60; % change units to minutes
        plot( 1:length(dimx_range) , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; % ? perhaps plot solid line only if number of thresholds exceeds 3?
        axis square
        xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
        xlabel(parx_lab) % set x axis label to match parameter name
        xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
        ylabel 'computation time (mins)'
        ylim([ min(o.m_tElapsed /60,[],'all')  max(o.m_tElapsed /60,[],'all')])
        ylim(ylim + [-range(ylim)/lms,range(ylim)/lms])
        title(subplot_title);
        
        % plot computation time (percentage)
        set(0, 'CurrentFigure', fig_ica_tElapsed_prop(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        eval([ 'temp2plot = o.m_tElapsed_prop(' indstr ');' ])
        plot( 1:length(dimx_range) , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; % ? perhaps plot solid line only if number of thresholds exceeds 3?
        axis square
        xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
        xlabel(parx_lab) % set x axis label to match parameter name
        xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
        ylabel 'computation time (% of raw)'
        ylim([ min(o.m_tElapsed_prop ,[],'all')  max(o.m_tElapsed_prop,[],'all')])
        ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
        title(subplot_title);
        
%         % plot identifiability (meanDiff)
%         set(0, 'CurrentFigure', fig_ica_identiDiff_meandiff(st)) % select this figure
%         subplot(nparr,nparc, getSubplotInd(r,c, nparc));
%         eval([ 'temp2plot = o.m_identiDiff_meanDiff(' indstr ');' ])
%         plot( 1:length(dimx_range) , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; % ? perhaps plot solid line only if number of thresholds exceeds 3?
%         axis square
%         xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
%         xlabel(parx_lab) % set x axis label to match parameter name
%         xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
%         ylabel 'mean identifiability score'
%         ylim([ -1  1])
%         title(subplot_title);
%         
%         % plot identifiability (propNumKnown)
%         set(0, 'CurrentFigure', fig_ica_identiDiff_propNumKnown(st)) % select this figure
%         subplot(nparr,nparc, getSubplotInd(r,c, nparc));
%         eval([ 'temp2plot = o.m_identiDiff_propNumKnown(' indstr ');' ])
%         plot( 1:length(dimx_range) , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; % ? perhaps plot solid line only if number of thresholds exceeds 3?
%         axis square
%         xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
%         xlabel(parx_lab) % set x axis label to match parameter name
%         xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
%         ylabel '% identifiable components'
%         ylim([ 0 100 ])
%         ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
%         title(subplot_title);
%         
%         % plot identifiability (propVarKnown)
%         set(0, 'CurrentFigure', fig_ica_identiDiff_propVarKnown(st)) % select this figure
%         subplot(nparr,nparc, getSubplotInd(r,c, nparc));
%         eval([ 'temp2plot = o.m_identiDiff_propVarKnown(' indstr ');' ])
%         plot( 1:length(dimx_range) , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; % ? perhaps plot solid line only if number of thresholds exceeds 3?
%         axis square
%         xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
%         xlabel(parx_lab) % set x axis label to match parameter name
%         xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
%         ylabel '% variance explained by identifiable components'
%         ylim([ 0 100 ])
%         ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
%         title(subplot_title);
        
        % plot class probabilities (% of components)
        set(0, 'CurrentFigure', fig_ica_classProb4_catNumProp(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        eval([ 'temp2plot = squeeze( o.m_classProb4_catNumProp(' indstr ',:) );' ])
        for k = 1:size(temp2plot,2)
            plot( 1:length(dimx_range) , temp2plot(:,k), '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; % ? perhaps plot solid line only if number of thresholds exceeds 3?
        end
        axis square
        xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
        xlabel(parx_lab) % set x axis label to match parameter name
        xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
        ylabel '% of components'
        ylim([0 100])
        ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
        if c == 1 && r == 1, legend(o.classProb4_classLabs, 'Location','northeast'); end
        title(subplot_title);
        
        % plot class probabilities (% of variance explained)
        set(0, 'CurrentFigure', fig_ica_classProb4_catVarProp(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        eval([ 'temp2plot = squeeze( o.m_classProb4_catVarProp(' indstr ',:) );' ])
        for k = 1:size(temp2plot,2)
            plot( 1:length(dimx_range) , temp2plot(:,k), '-o', 'LineWidth', 1.5, 'MarkerSize', MarkerSize ); hold on; % ? perhaps plot solid line only if number of thresholds exceeds 3?
        end
        axis square
        xticks(1:length(dimx_range)); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
        xlabel(parx_lab) % set x axis label to match parameter name
        xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
        ylabel '% of variance explained'
        ylim([0 100])
        ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
        if c == 1 && r == 1, legend(o.classProb4_classLabs, 'Location','northeast'); end
        title(subplot_title);
    
    
        drawnow
    end
    
    %% end loop through parameters plotted on seperate figures
        end % end row loop
    end % end col loop

    %% end loop through stages


end    % end stage loop

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
    figs = who('fig_*'); % get all figure handles
    for f = 1:length(figs)
        eval([ 'savefig(' figs{f}  ', [''fig_' cfg.saveAs '_'  figs{f}(5:end) '''])' ]) % save fig
    end
    
    fprintf('... saved figures\n')
end

%% Subfunctions
function ind = getSubplotInd(r,c,numc)
    ind = (r-1)*numc + c;
end


end % END MAIN FUNCTION

