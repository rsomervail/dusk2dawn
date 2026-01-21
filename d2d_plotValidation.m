% 
% 
%       plot validation metrics from up to 3 parameter dimensions at a time (for one dataset or a group)
%
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net 
%%  
function d2d_plotValidation(EEG, cfg)

%% handle multiple datasets
nfiles = length(EEG);

% check group-level merging has been done, if not then do it now
if nfiles > 1
    EEG = d2d_group_validateMerge(EEG);

    fprintf('\n%s: plotting effects of ASR on various validation metrics for the following datasets:\n',mfilename)
    disp({EEG(:).setname}')
%     disp(EEG(1).etc.dusk2dawn.group.datasets') 
else
    fprintf('\n%s: plotting effects of ASR on various validation metrics for dataset:\n%s\n',mfilename, EEG.setname)
end

%% defaults
% general
if ~isfield(cfg,'saveFig'),     cfg.saveFig   = false;   end
if ~isfield(cfg,'maxFig'),      cfg.maxFig   = true;   end
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
LineWidth  = 1.5;
AreaAlpha  = 0.5; 
SingleSubOpacity = 0.4;

cols_classProb4 = [ 0,0,1; 1,0.5,0; 0,0,0; 1,0,0  ]; % brain=blue, bioArtifact=orange, recArtifact=black, unknown=red
col_orange = [1,0.4,0];
cols_freqs = [0,0,1;1,0,0;0,1,0;0,0,0.172413793103448;1,0.103448275862069,0.724137931034483;1,0.827586206896552,0;0,0.344827586206897,0;0.517241379310345,0.517241379310345,1;0.620689655172414,0.310344827586207,0.275862068965517;0,1,0.758620689655172];

%% get dataset info
cfgd = EEG(1).etc.dusk2dawn.cfg;

%% get parameters to plot on each axis
pars  = cfgd.pars;
npars = length(pars.labels);

% if no varied parameters then create virtual parameter and set this to plot on X axis
if npars == 0
    pars.labels = {'ASR'};
    pars.values{1} = {'cleaned'};
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
parx_val = [ 'raw'; parx_val  ]; % add raw data value
parx_val = regexprep(parx_val, '\s*',' '); % trim excess spaces#

% get 2nd and 3rd cleaning parameter dimensions (for subplots)
if length( cfg.dims2plot ) >= 2
    dimc = cfg.dims2plot(2); 
    parc_lab = pars.labels{dimc};
    parc_lab = strrep(parc_lab,'_','-'); % replace underscores for plotting
    parc_val = d2d_getParList(pars.values{dimc}); % extract varied parameter values
    parc_val = cellfun(@num2str, parc_val, 'UniformOutput', false);
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
if ~cfgd.splitByStage
    stageCodes = {''};
    nstages = 1;
    cfg.stages = 1;
else
    stageCodes = cfgd.stageCodes;
    stageCodes = stageCodes(cfg.stages);
    nstages = length(stageCodes);
end

for st = 1:nstages

    %% get validation metrics
    if nfiles == 1
        origFile = strrep(EEG.etc.dusk2dawn.cfg.origFile,'.set','');
        v = EEG.etc.dusk2dawn.valid_merged(st);
    else
        v = EEG(1).etc.dusk2dawn.group.validG(cfg.stages(st)); % every dataset has the same group-level validation structure
    end

    %% Initialise all figures
    if nfiles == 1
        if cfgd.splitByStage
            figprefix = [ 'D2D - ' origFile '_' stageCodes{st} ': ' ];
        else
            figprefix = [ 'D2D - ' origFile ': ' ];
        end
    else
        if cfgd.splitByStage
            figprefix = [  'D2D - group_' stageCodes{st} ': ' ];
        else
            figprefix = [  'D2D - group: '];
        end
    end

    % cleanData
    if cfg.plot_asr
        fig_asr_tElapsed_abs(st)  = figure('name',[ figprefix  'ASR computation time?'], 'numbertitle','off');
        fig_cleanData(st)         = figure('name',[ figprefix 'how much variance did ASR remove?'], 'numbertitle','off');
        fig_calibData(st)         = figure('name',[ figprefix  'how much data did ASR have for calibration?'], 'numbertitle','off');
    end
    % fft
    if isfield(v,'fft') && cfg.plot_fft
        fig_fft_powerByBand_abs(st)  = figure('name',[ figprefix  'FFT amplitude by band (absolute)'], 'numbertitle','off');
        fig_fft_powerByBand_prop(st) = figure('name',[ figprefix  'FFT amplitude by band (%)'], 'numbertitle','off');
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
        fig_ica_classProb4_catNumProp(st)      =  figure('name',[ figprefix  'ICA % ICs of each category'], 'numbertitle','off'); %#ok<*AGROW> 
        fig_ica_classProb4_catVarProp(st)      =  figure('name',[ figprefix  'ICA % variance explained by each category'], 'numbertitle','off');
    end
    
    %% Loop through cleaning parameters which are plotted on seperate figures
    % get range to loop through for the row and column dimensions (x axis is handled below)
    if isnan(dimr) 
        dimr_range = 1; % if not even plotting this dimension, just loop through 1
    elseif strcmp(cfg.plotRange{dimr},':')
        dimr_range = 1:nparr; % use entire range of values
    else
        eval([ 'dimr_range = ' cfg.plotRange{dimr} ';' ]) %#ok<*EVLEQ> % use subset of values  
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
    nparx = length(dimx_range);

    for r = 1:nparr
        for c = 1:nparc
            
    %% construct a plot-index string for accessing validation matrices specific for this iteration of the row/column loop
    indstr = cfg.plotRange; % get plot range for x axis (row and column ranges will be overwritten by the values for this iteration)
    if ~isnan(dimc), indstr{dimc} = num2str(dimc_range(c)+1); end % after these two steps, should be left with one index corresponding to the x axis range
    if ~isnan(dimr), indstr{dimr} = num2str(dimr_range(r)+1); end % +1 because first element of each dim is always original data
    if nfiles == 1
        indstr =       strjoin(indstr,',');  % make into string, first colon refers to dimension of merged datasets
    else
        indstr = [':,' strjoin(indstr,',')]; % make into string, first colon refers to dimension of merged datasets
    end
    
    %% construct title for each subplot
    switch length( cfg.dims2plot )
        case 1
            subplot_title = [];
        case 2
            subplot_title = [parc_lab ': ' parc_val{dimc_range(c)}  ];
        case 3
            subplot_title = [parc_lab ': ' parc_val{dimc_range(c)} '   ' parr_lab ': ' parr_val{dimr_range(r)}   ];
    end
    
    %% PLOT - ASR cleaning metrics
    if isfield(v,'asr') && cfg.plot_asr
        o = v.asr; % get subfield
        
        % computation time (ASR) (absolute)
        set(0, 'CurrentFigure', fig_asr_tElapsed_abs(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        met2plot = 'm_tElapsed';
        eval([ 'temp2plot = o.(met2plot)(' indstr ');' ])
        temp2plot = double(temp2plot);
        temp2plot = temp2plot ./ 60; % change units to minutes
        if nfiles == 1
            plot( dimx_range , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on;
        else
            temp2plot = squeeze(temp2plot);
            [lineh, areah] = boundedline( dimx_range, mean(temp2plot), std(temp2plot) ,'-or'  ); hold on
        end
        ylabel 'computation time (mins)'
        setYlim( 'abs' , 1/60 ) % plot type, multiplier for data
            editPlot; % set various plotting parameters
            reorderPlots(gca); % bring all lines to foreground and send patches to background
    
        % percentage of timepoints altered by ASR
        set(0, 'CurrentFigure', fig_cleanData(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        met2plot = 'm_propCleaned';
        yyaxis left
        eval([ 'temp2plot = o.(met2plot)(' indstr ');' ]) 
        temp2plot = double(temp2plot);
        if nfiles == 1
            plot( dimx_range , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; 
        else
            temp2plot = squeeze(temp2plot);
            [lineh, areah] = boundedline( dimx_range, mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
        end
        ylabel '% timepoints altered'
        setYlim('%',1) % plot type, multiplier
            editPlot; % set various plotting parameters
            reorderPlots(gca); % bring all lines to foreground and send patches to background
        % percentage of variance removed by ASR
        yyaxis right
        met2plot = 'm_propVarRemoved';
%         met2plot = 'm_propGFPRemoved';
        eval([ 'temp2plot = o.(met2plot)(' indstr ');' ]) 
        temp2plot = double(temp2plot);
        if nfiles == 1
            plot( dimx_range , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize, 'Color',col_orange ); hold on; 
        else
            temp2plot = squeeze(temp2plot);
            [lineh, areah] = boundedline( dimx_range, mean(temp2plot), std(temp2plot) ,'-o','cmap',col_orange  ); hold on
        end
        ylabel '% variance removed'
        setYlim('%',1) % plot type, multiplier
            editPlot; % set various plotting parameters
            reorderPlots(gca); % bring all lines to foreground and send patches to background
        
        % how sufficient was calibration data?
        % minimum calibration times
        set(0, 'CurrentFigure', fig_calibData(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        yyaxis left
        met2plot = 'm_calibDataLen_min';
        eval([ 'temp2plot = o.(met2plot)(' indstr ');' ]) 
        temp2plot = double(temp2plot);
        if nfiles == 1
            plot( dimx_range , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on;
        else
            temp2plot = squeeze(temp2plot);
            [lineh, areah] = boundedline( dimx_range, mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
            for sb = 1:size(temp2plot,1), plot( dimx_range, temp2plot(sb,:), '-', 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
        end
        ylabel 'minimum calibration time (s)'
        setYlim('abs',1)
        plot(xlim,[60,60],'r--'); % 60s threshold dashed red line
            editPlot; % set various plotting parameters
            reorderPlots(gca); % bring all lines to foreground and send patches to background
        % percentage of data chunks with sufficient calibration
        yyaxis right
        met2plot = 'm_calibDataLen_prop60s';
        eval([ 'temp2plot = o.(met2plot)(' indstr ');' ])
        temp2plot = double(temp2plot);
        if nfiles == 1
            plot( dimx_range , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; 
        else
            temp2plot = squeeze(temp2plot);
            [lineh, areah] = boundedline( dimx_range, mean(temp2plot), std(temp2plot) ,'-o', 'cmap',[1,0.4,0]  ); hold on
        end
        ylabel '% of cleaning chunks with sufficient calibration data (>=60s)'
        setYlim('%',1)
            editPlot; % set various plotting parameters
            reorderPlots(gca); % bring all lines to foreground and send patches to background
    
        drawnow
    end 
    %% PLOT - frequency decomposition metrics
    if isfield(v,'fft') && cfg.plot_fft
        o = v.fft; % get subfield
    
        % plot frequency power by band (absolute)
        set(0, 'CurrentFigure', fig_fft_powerByBand_abs(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        met2plot = 'm_binAmp_avgChan';
        eval([ 'temp2plot = squeeze( o.(met2plot)(' indstr ',:) );' ])
        temp2plot = double(temp2plot);
        if nfiles == 1
            for k = 1:size(temp2plot,2)
                lineh(k) = plot( dimx_range , temp2plot(:,k), '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on;
            end
        else
            temp2plot = squeeze(temp2plot);
            clear lineh areah
            for k = 1:size(temp2plot,3)
                [lineh(k), areah(k)] = boundedline( dimx_range, mean(temp2plot(:,:,k)), std(temp2plot(:,:,k)) ,'-o', 'cmap', cols_freqs(k,:)  ); hold on
            end
        end
        ylabel 'amplitude (normalisd by frequency)'
        setYlim('abs',1)
        if c == 1 && r == 1, legend(lineh,o.binFreqsLabels, 'Location','northeast'); end
            editPlot; % set various plotting parameters
            reorderPlots(gca); % bring all lines to foreground and send patches to background
        
        
        % plot frequency power by band (percentage)
        set(0, 'CurrentFigure', fig_fft_powerByBand_prop(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        met2plot = 'm_binAmp_avgChan_prop';
        eval([ 'temp2plot = squeeze( o.(met2plot)(' indstr ',:) );' ])
        temp2plot = double(temp2plot);
        if nfiles == 1
            for k = 1:size(temp2plot,2)
                lineh(k) = plot( dimx_range , temp2plot(:,k), '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; 
            end
        else
            temp2plot = squeeze(temp2plot);
            clear lineh areah
            for k = 1:size(temp2plot,3)
                [lineh(k), areah(k)] = boundedline( dimx_range, mean(temp2plot(:,:,k)), std(temp2plot(:,:,k)) ,'-o', 'cmap', cols_freqs(k,:)  ); hold on
            end
        end
        ylabel 'amplitude (%)'
        setYlim('abs',1) % abs here because makes a neater plot for power reduction
        if c == 1 && r == 1, legend(lineh,o.binFreqsLabels, 'Location','northeast'); end
            editPlot; % set various plotting parameters
            reorderPlots(gca); % bring all lines to foreground and send patches to background
       
        drawnow
    end
    %% PLOT - slow-wave metrics
    if isfield(v,'sw') && cfg.plot_sw
        o = v.sw; % get subfield
    
        % plot n slow-waves
        set(0, 'CurrentFigure', fig_sw_nwaves(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        met2plot = 'm_nwaves';
        eval([ 'temp2plot = o.(met2plot)(' indstr ');' ])
        temp2plot = double(temp2plot);
        if nfiles == 1
            plot( dimx_range , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; 
        else
            temp2plot = squeeze(temp2plot);
            [lineh, areah] = boundedline( dimx_range, mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
            for sb = 1:size(temp2plot,1), plot( dimx_range, temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
        end
        ylabel 'number of slow-waves'
        setYlim('abs',1);
            editPlot; % set various plotting parameters
            reorderPlots(gca); % bring all lines to foreground and send patches to background
        
        % plot slow-wave amplitude (absolute)
        set(0, 'CurrentFigure', fig_sw_amp_mean_abs(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        met2plot = 'm_amp_mean';
        eval([ 'temp2plot = o.(met2plot)(' indstr ');' ])
        temp2plot = double(temp2plot);
        if nfiles == 1
            plot( dimx_range , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; 
        else
            temp2plot = squeeze(temp2plot);
            [lineh, areah] = boundedline( dimx_range, mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
            for sb = 1:size(temp2plot,1), plot( dimx_range, temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
        end
        ylabel 'slow-wave amplitude (µV)'
        setYlim('abs',1);
            editPlot; % set various plotting parameters
            reorderPlots(gca); % bring all lines to foreground and send patches to background
        
        % plot slow-wave amplitude (percentage)
        set(0, 'CurrentFigure', fig_sw_amp_mean_prop(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        met2plot = 'm_amp_mean_prop';
        eval([ 'temp2plot = o.(met2plot)(' indstr ');' ])
        temp2plot = double(temp2plot);
        if nfiles == 1
            plot( dimx_range , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; 
        else
            temp2plot = squeeze(temp2plot);
            [lineh, areah] = boundedline( dimx_range, mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
            for sb = 1:size(temp2plot,1), plot( dimx_range, temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
        end
        ylabel 'slow-wave amplitude (% of raw)'
        setYlim('%',1);
            editPlot; % set various plotting parameters
            reorderPlots(gca); % bring all lines to foreground and send patches to background
        
        % plot slow-wave t-value
        set(0, 'CurrentFigure', fig_sw_tval_median(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        met2plot = 'm_tval_median';
        eval([ 'temp2plot = o.(met2plot)(' indstr ');' ])
        temp2plot = double(temp2plot);
        if nfiles == 1
            plot( dimx_range , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; 
        else
            temp2plot = squeeze(temp2plot);
            [lineh, areah] = boundedline( dimx_range, mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
            for sb = 1:size(temp2plot,1), plot( dimx_range, temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
        end
        ylabel 'slow-wave consistency (t value)'
        setYlim('abs',1);
            editPlot; % set various plotting parameters
            reorderPlots(gca); % bring all lines to foreground and send patches to background
    
        drawnow
    end
    
    %% PLOT - ica metrics
    if isfield(v,'ica') && cfg.plot_ica
        o = v.ica; % get subfield
    
        % plot computation time (absolute)
        set(0, 'CurrentFigure', fig_ica_tElapsed_abs(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        met2plot = 'm_tElapsed';
        eval([ 'temp2plot = o.(met2plot)(' indstr ');' ])
        temp2plot = double(temp2plot);
        temp2plot = temp2plot ./ 60; % change units to minutes
        if nfiles == 1
            plot( dimx_range , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; 
        else
            temp2plot = squeeze(temp2plot);
            [lineh, areah] = boundedline( dimx_range, mean(temp2plot), std(temp2plot) ,'-or'  ); hold on
            for sb = 1:size(temp2plot,1), plot( dimx_range, temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
        end
        ylabel 'computation time (mins)'
        setYlim('abs',1/60);
            editPlot; % set various plotting parameters
            reorderPlots(gca); % bring all lines to foreground and send patches to background
        
        % plot computation time (percentage)
        set(0, 'CurrentFigure', fig_ica_tElapsed_prop(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        met2plot = 'm_tElapsed_prop';
        eval([ 'temp2plot = o.(met2plot)(' indstr ');' ])
        temp2plot = double(temp2plot);
        if nfiles == 1
            plot( dimx_range , temp2plot, '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize ); hold on; 
        else
            temp2plot = squeeze(temp2plot);
            [lineh, areah] = boundedline( dimx_range, mean(temp2plot), std(temp2plot) ,'-o'  ); hold on
            for sb = 1:size(temp2plot,1), plot( dimx_range, temp2plot(sb,:), 'Color', [0.2, 0.2, 0.2, SingleSubOpacity]  ); end % plot single subjects
        end
        ylabel 'computation time (% of raw)'
        setYlim('abs',1) % more zoomed in for subtle differences in % computation time
            editPlot; % set various plotting parameters
            reorderPlots(gca); % bring all lines to foreground and send patches to background
        
        % don't have identifiability plots here yet (low priority)
        
        % plot class probabilities (% of components)
        set(0, 'CurrentFigure', fig_ica_classProb4_catNumProp(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        met2plot = 'm_classProb4_catNumProp';
        eval([ 'temp2plot = squeeze( o.(met2plot)(' indstr ',:) );' ])
        temp2plot = double(temp2plot);
        clear lineh areah
        if nfiles == 1
            for k = 1:size(temp2plot,2)
                lineh(k) = plot( dimx_range , temp2plot(:,k), '-o', 'LineWidth', 1.5, 'MarkerSize',MarkerSize, 'Color', cols_classProb4(k,:)  ); hold on; % ? perhaps plot solid line only if number of thresholds exceeds 3?
            end
        else
            temp2plot = squeeze(temp2plot);
            for k = 1:size(temp2plot,3)
                [lineh(k), areah(k)] = boundedline( dimx_range, mean(temp2plot(:,:,k)), std(temp2plot(:,:,k)) ,'-o', 'cmap', cols_classProb4(k,:)  ); hold on
            end
        end
        ylabel '% of components'
        setYlim('%',1)
        if c == 1 && r == 1, legend(lineh,o.classProb4_classLabs, 'Location','northeast'); end
            editPlot; % set various plotting parameters
            reorderPlots(gca); % bring all lines to foreground and send patches to background
        
        % plot class probabilities (% of variance explained)
        set(0, 'CurrentFigure', fig_ica_classProb4_catVarProp(st)) % select this figure
        subplot(nparr,nparc, getSubplotInd(r,c, nparc));
        met2plot = 'm_classProb4_catVarProp';
        eval([ 'temp2plot = squeeze( o.(met2plot)(' indstr ',:) );' ])
        temp2plot = double(temp2plot);
        clear lineh areah
        if nfiles == 1
            for k = 1:size(temp2plot,2)
                lineh(k) = plot( dimx_range , temp2plot(:,k), '-o', 'LineWidth', 1.5, 'MarkerSize', MarkerSize, 'Color', cols_classProb4(k,:) ); hold on; % ? perhaps plot solid line only if number of thresholds exceeds 3?
            end
        else
            temp2plot = squeeze(temp2plot);
            for k = 1:size(temp2plot,3)
                [lineh(k), areah(k)] = boundedline( dimx_range, mean(temp2plot(:,:,k)), std(temp2plot(:,:,k)) ,'-o', 'cmap', cols_classProb4(k,:)  ); hold on
            end
        end
        axis square
        ylabel '% of variance explained'
        setYlim('%',1)
        if c == 1 && r == 1, legend(lineh,o.classProb4_classLabs, 'Location','northeast'); end
            editPlot; % set various plotting parameters
            reorderPlots(gca); % bring all lines to foreground and send patches to background
    
        drawnow
    end
    
        end % end row loop
    end % end col loop
end    % end stage loop

%% Maximise all figures produced
if cfg.maxFig
    figs = findall(0, 'type', 'Figure');
    figs(~startsWith({figs.Name},'D2D - ')) = [];
    for f = 1:length(figs)
        fig = figs(f);
        fig.Units = 'normalized';
        fig.Position = [0 0.0370 1 0.8917];
    end
end

%% Save figures
if cfg.saveFig
    fprintf('saving figures ...\n')

    % get savedir
    cfg.savePath = EEG(1).etc.dusk2dawn.cfg.savePath;

    % loop through figure handles and save
    figvars = whos('fig_*');
    fighandles = nan(size(figvars,1) * nstages,1);
    fignames   = cell(size(figvars,1) * nstages,1);
    count = 1;
    for st = 1:nstages
        stage = stageCodes{st};
        stage = strrep(stage,'[','');
        stage = strrep(stage,']','');
        stage = [stage '_'];
        for f = 1:length(figvars)
            fignames{count} = strrep(figvars(f).name,'fig_',['fig_' stage]);
            eval([  'fighandles(count) = ' figvars(f).name '(st);'])
            count = count + 1;
        end
    end
    for f = 1:length(fignames)
        savefig( fighandles(f) , [cfg.savePath  filesep  fignames{f} ]); % save fig
    end
    
    fprintf('... saved figures\n')
end

%% Subfunctions
function setYlim(mode, mult)
    switch mode 
        case '%'
            ylim([0 100]) % already expanded in editPlot function
%             ylim([-10 110])
        case 'abs'
            if range(o.(met2plot),'all') ~= 0
                ylim([ min(o.(met2plot) * mult,[],'all')  max(o.(met2plot) * mult,[],'all')])
            end
    end
end
function editPlot
%     ax = gca;
    % group vs single-dataset specific
    if nfiles > 1 % group-level
        % set linewidths and such
        for lh = 1:length(lineh)
            lineh(lh).LineWidth         = LineWidth;
            lineh(lh).MarkerSize        = MarkerSize;
            lineh(lh).MarkerFaceColor   = lineh(lh).Color;
            areah(lh).FaceAlpha         = AreaAlpha;
        end
    else % single-dataset
        % ? might be nothing that needs to go here
    end
    % general
    axis square
    title(subplot_title);
    xticks(dimx_range); xticklabels(parx_val(dimx_range)) % set x axis ticks to match parameter values
    xlabel(parx_lab) % set x axis label to match parameter name
    xlim(xlim + [-range(xlim)/lms,range(xlim)/lms]) % expand x-axis a bit
    ylim(ylim + [-range(ylim)/lms,range(ylim)/lms]) % expand y-axis a bit
end
function reorderPlots(ax)
    % assuming patch is always last
    lines = []; theRest = []; grouplines = [];
    for j = 1 :length(ax.Children)
        if strcmp( ax.Children(j).Type, 'line' )
            if ~strcmp( ax.Children(j).Marker , 'none' ) % average line
                grouplines = [grouplines, j];
            else % single subject line
                lines = [lines, j];
            end
        else % everything else (patches for error bars) goes to the background
            theRest = [theRest, j];
        end
    end
    % change plot order
    ax.Children = ax.Children( [ grouplines, lines, theRest ]  );   
end

function ind = getSubplotInd(r,c,numc)
    ind = (r-1)*numc + c;
end


end % END MAIN FUNCTION

