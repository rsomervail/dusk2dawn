%
%
%       - wrapper function for ASR calibration and cleaning
%       - takes one EEGlab dataset
%       - can calibrate & clean the whole dataset together, or slide through
%          the dataset in shorter chunks to calibrate according to chunk-specific thresholds
% 
%
% 
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net
% History:
% 19/01/2023 ver 1.0.0 Created
% 
%%  
function [EEG_out, cfg] = d2d_cleanData(EEG, cfg)

% defaults
if ~isfield(cfg,'splitBySlidingWindow'), cfg.splitBySlidingWindow   = false; end
if ~isfield(cfg,'chunk_len'),     cfg.chunk_len                 = inf;   end
if ~isfield(cfg,'chunk_overlap'), cfg.chunk_overlap             = 1;   end
if ~isfield(cfg,'checkCalibManual'),   cfg.checkCalibManual     = false; end
if ~isfield(cfg,'checkCalibManual_plotStyle'),   cfg.checkCalibManual_plotStyle = 'butterfly'; end
if ~isfield(cfg,'checkCalibManual_channels'),  cfg.checkCalibManual_channels = []; end
if ~isfield(cfg,'ref_maxbadchannels'),   cfg.ref_maxbadchannels   = false; end
if ~isfield(cfg,'ref_tolerances'),   cfg.ref_tolerances   = false; end
if ~isfield(cfg,'ref_wndlen'),   cfg.ref_wndlen   = []; end
if ~isfield(cfg,'asr_useRiemmanian'),   cfg.asr_useRiemmanian   = false; end
if ~isfield(cfg,'asr_maxdims'),   cfg.asr_maxdims   = 2/3; end
if ~isfield(cfg,'asr_useGPU'),   cfg.asr_useGPU   = false; end
if ~isfield(cfg,'asr_MaxMem')   
    try 
        cfg.asr_MaxMem   = d2d_recMaxMem(cfg.asr_useGPU); 
    catch
        cfg.asr_MaxMem = [];
    end
end


%% load dataset from disk if no variable was passed in from memory  
if isempty(EEG)
    EEG = d2d_loadDataset(cfg);
end

fprintf('d2d_cleanData: cleaning dataset "%s" ... \n', EEG.setname) 

% remove info from raw data
EEG.etc.dusk2dawn = rmfield( EEG.etc.dusk2dawn, 'valid' );

%% if cleaning in chunks (i.e. sliding a window through data and applying calibration separately to each segment)
% turn this off if chunk length is infinite
if ~isfinite(cfg.chunk_len), cfg.splitBySlidingWindow  = false; end
if cfg.splitBySlidingWindow 
    %% misc stuff
    % check options make sense
    if cfg.chunk_len < 1
        error 'd2d_cleanData: minimum reasonable chunk length is at least 1 minute (in order to have enough calibration data)'; 
    end
    if      cfg.chunk_overlap >= cfg.chunk_len
        error 'd2d_cleanData: overlap of cleaning chunks cannot be longer than the chunks themselves';
    end
    
    % currently doesn't support manual calibration checking 
    cfg.checkCalibManual = false;
    
    %% define chunks   
    nsamps   = length(EEG.times); % get number of samples    
    chunklen = cfg.chunk_len     * 60 * EEG.srate; % calculate chunk length in timepoints
    chunklap = cfg.chunk_overlap * 60 * EEG.srate; % get chunk overlap in timepoints
    chunks(1,1) = 1;
    chunks(1,2) = 1 + chunklen;
    while chunks(end,2) < nsamps 
        chunks(end+1,1) = chunks(end,2) + 1 - chunklap; %#ok<AGROW>
        chunks(end  ,2) = chunks(end,1) + chunklen;
    end
    chunks(end,end) = nsamps; % last chunk will have to go to the end
    nchunks = size(chunks,1);
    
    % outputs
    cfgout.chunks  = chunks;
    cfgout.nchunks = nchunks;
    
    % prepare chunk variable
    EEG_chunk = EEG; EEG_chunk.data = [];
    EEG_chunk.pnts = chunklen + 1;
    
    %% create output structure for data
    EEG_out = EEG;
    EEG_out.data = nan(size(EEG_out.data));
    
    %% loop through chunks and apply ASR cleaning 
    fprintf('running ASR cleaning ...\n')
    tin_asr = tic;
    for ch = 1 : nchunks
        
        %% get chunk
        EEG_chunk.data = EEG.data( :, chunks(ch,1):chunks(ch,2) );

        %% get calibration data        
        evalc([ ' [ref_section, calibMask] = clean_windows(EEG, cfg.ref_maxbadchannels, cfg.ref_tolerances, cfg.ref_wndlen);   ' ]);
        
        %% check calibration data are an adequate length
        calibLen = size(ref_section.times,2) / ref_section.srate; % compute calibration data length in seconds
        fprintf('d2d_cleanData: using %.2f mins of calibration data ...\n', calibLen/60)
        if      calibLen < 60
            warning('d2d_cleanData: chosen calibration data are shorter than a minute, this is not ideal -> consider increasing chunk length')
        elseif  calibLen < 15
            error('Error: chosen calibration data are shorter than 15 seconds, this is not adequate -> increase chunk length')
        end
        cfgout.calibDataLen(ch)    = calibLen; 
        cfgout.calibDataMask{ch}   = calibMask;

        %% run ASR cleaning 
        EEG_chunk_cleaned = clean_asr( EEG_chunk, cfg.asr_cutoff, [],[],cfg.asr_maxdims,     ...
                   ref_section,[],[], cfg.asr_useGPU , cfg.asr_useRiemmanian, cfg.asr_MaxMem ); 
        
        %% insert chunk into output data structure
        % take average of this cleaned chunk and any segment of overlapping data that may have already been processed
        EEG_out.data(:, chunks(ch,1):chunks(ch,2) ) = mean( cat( 3, EEG_out.data(:, chunks(ch,1):chunks(ch,2) ), EEG_chunk_cleaned.data ), 3 ,'omitnan'); 
        
        fprintf('chunk %d/%d (%.1f%%) completed ...\n\n', ch, nchunks, 100*ch/nchunks )
    end % chunk loop
    cfgout.tElapsed = toc(tin_asr);
    fprintf( '... cleaning completed in %.2f mins\n\n', cfgout.tElapsed/60) 
    
    %% compute general ASR stats 
    % find altered timepoints and compute proportion 
    cfgout.pointsCleaned = any( EEG_out.data ~= EEG.data );  
    cfgout.propCleaned =  100 * sum(cfgout.pointsCleaned) / EEG.pnts; 

    % compute variance removed (% of GFP summed across time) 
    preVar = sum(GFP( EEG.data-mean(EEG.data) )); % avg ref before computing GFP
    posVar = sum(GFP( EEG_out.data-mean(EEG_out.data)           ));
    cfgout.propGFPRemoved = 100 * ( preVar - posVar ) / preVar; clearvars preVar posVar

    % compute variance removed (% of summed absolute signal across channels and time)  ? what's the best metric to use 
    preVar = sum(abs( EEG.data )      ,'all'); 
    posVar = sum(abs( EEG_out.data )  ,'all');
    cfgout.propVarRemoved = 100 * ( preVar - posVar ) / preVar; clearvars preVar posVar
    
%% if cleaning whole dataset in one chunk
else  
    
    %% find initial calibration data  
    if isfield(cfg,'calibrationData') % check if this has been provided already by the user
        fprintf('d2d_cleanData: using user specified calibration data\n')
        ref_section = cfg.calibrationData;
    else % if no calibration data is provided by the user
        fprintf('d2d_cleanData: finding data to use for calibration ...\n')
        evalc([ ' [ref_section, calibMask] = clean_windows(EEG, cfg.ref_maxbadchannels, cfg.ref_tolerances, cfg.ref_wndlen);   ' ]);

        % if manual check is active, edit the calibration selection
        if cfg.checkCalibManual
            fprintf('d2d_cleanData: plotting the automatically chosen calibration data -> check that the selected data are relatively clean\n (clear segments) and mark artifacts segments (in red)\n')
            
            % choose channels
            if isempty(cfg.checkCalibManual_channels), cfg.checkCalibManual_channels = {EEG.chanlocs.labels}; end
            if iscell(cfg.checkCalibManual_channels)
                allChans = {EEG.chanlocs.labels};
                [~, cfg.checkCalibManual_channels] = ismember(cfg.checkCalibManual_channels,allChans); clear allChans
            end
            EEG2plot            = EEG; 
            EEG2plot.chanlocs   = EEG.chanlocs(cfg.checkCalibManual_channels); 
            EEG2plot.data       = EEG.data(cfg.checkCalibManual_channels,:);
            EEG2plot.nbchan     = length(cfg.checkCalibManual_channels);
            
            % convert to fieldtrip and plot
            EEG2plot_ft = eeglab2fieldtrip( EEG2plot, 'raw', 'none');
            ctemp = [];
            ctemp.artfctdef.visual.artifact = logical2indices(~calibMask); % load in the automatically found calibration data
            ctemp.continuous = 'yes';
            switch cfg.checkCalibManual_plotStyle   % ! add to documentation
                case 'butterfly'     
                    ctemp.viewmode =  'butterfly'; ctemp.ylim = [-200 200];
                case 'vertical'
                    ctemp.viewmode =  'vertical'; ctemp.ylim  = [-(480.1*EEG2plot.nbchan^-0.8802) (480.1*EEG2plot.nbchan^-0.8802)]; 
            end
            ctemp.blocksize  = 30; 
            ctemp.continuous   = 'yes';
                cfglay.elec = EEG2plot_ft.elec;
            ctemp.layout = ft_prepare_layout( cfglay ); clear cfglay
            ctemp = ft_databrowser(ctemp, EEG2plot_ft); % plot data and allow user to choose calibration data
            
            % get chosen calibration data
            calibMask = ~indices2logical( ctemp.artfctdef.visual.artifact, length(calibMask) ); % take the points which are NOT marked as artifacts
            ref_section = pop_select(EEG, 'point', logical2indices(calibMask)); 
        end
    end % find calibration data
    
    %% generate logical mask of calibration data (if not generated in previous section) (for later validation & analysis)
    if ~exist('calibMask','var')
        calibMask = ref_section.etc.clean_sample_mask;
    end
    
    %% check calibration data are adequate length
    calibLen = size(ref_section.times,2) / ref_section.srate; % compute calibration data length in seconds
    fprintf('... using %.2f mins of data for calibration\n', calibLen/60)
    if      calibLen < 60
        warning('chosen calibration data are shorter than a minute, this is not ideal -> manually check to see if more relatively clean data can be included')
    elseif  calibLen < 15
        error('Error: chosen calibration data are shorter than 15 seconds, this is not adequate -> manually check to see if more relatively clean data can be included')
    end
    cfgout.calibDataLen    = calibLen; 
    cfgout.calibDataMask   = calibMask;

    %% run ASR cleaning 
    fprintf('running ASR cleaning ...\n')
    tin_asr = tic;
        EEG_out = clean_asr( EEG, cfg.asr_cutoff, [],[],cfg.asr_maxdims,     ...
                       ref_section,[],[], cfg.asr_useGPU , cfg.asr_useRiemmanian, cfg.asr_MaxMem ); 
    cfgout.tElapsed = toc(tin_asr);
    fprintf( '... cleaning completed in %.2f mins\n\n', cfgout.tElapsed/60 ) 
    
    %% compute general ASR stats
    % find altered timepoints and compute proportion 
    cfgout.pointsCleaned = any( EEG_out.data ~= EEG.data );  % ? this is very light rel to real data so storing is fine and worthwhile
    cfgout.propCleaned =  100 * sum(cfgout.pointsCleaned) / EEG.pnts; 

    % compute variance removed (% of GFP summed across time) 
    preVar = sum(GFP( EEG.data-mean(EEG.data) )); % avg ref before computing GFP
    posVar = sum(GFP( EEG_out.data-mean(EEG_out.data)           ));
    cfgout.propGFPRemoved = 100 * ( preVar - posVar ) / preVar; clearvars preVar posVar
    
    % compute variance removed (% of summed absolute signal across channels and time)  ? what's the best metric to use 
    preVar = sum(abs( EEG.data )      ,'all'); 
    posVar = sum(abs( EEG_out.data )  ,'all');
    cfgout.propVarRemoved = 100 * ( preVar - posVar ) / preVar; clearvars preVar posVar
    
end % end if cleaning in chunks

%% outputs
% store ouputs and settings in EEG data
cfgout.raw = false;
cfgout.pars = cfg.pars;
EEG_out.etc.dusk2dawn.asr = cfgout;

%% save dataset in specified path
if isfield(cfg, 'savePath')
    d2d_saveDataset(EEG_out,cfg);
end


end