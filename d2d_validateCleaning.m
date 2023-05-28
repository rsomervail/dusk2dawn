% 
%   - computes validation metrics from a dataset processed by Dusk2Dawn
%       - data can be raw or ASR-cleaned (with D2D)
%
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net 
%%  
function [EEG, valid] = d2d_validateCleaning(EEG, cfg)

% defaults - FFT
if ~isfield(cfg, 'fft'), cfg.fft = struct; end
if ~isfield(cfg.fft, 'run'), cfg.fft.run = true; end
if ~isfield(cfg.fft, 'binFreqs'), cfg.fft.binFreqs = [ 1,4;  4,8;  8,12; 12,16;  18,30;   30,45;  45,100  ]; end
if ~isfield(cfg.fft, 'binFreqsLabels'), cfg.fft.binFreqsLabels = { 'delta'; 'theta'; 'alpha'; 'sigma'; 'beta'; 'gammaLow'; 'gammaHigh'  }; end
if ~isfield(cfg.fft, 'storeWholeFFT'), cfg.fft.storeWholeFFT = true; end

% defaults - SW
if ~isfield(cfg, 'sw'), cfg.sw = struct; end
if ~isfield(cfg.sw, 'run'), cfg.sw.run = true; end 
if ~isfield(cfg.sw, 'event'), cfg.sw.event = []; end  % later this is checked to decide whether to auto-find slow-wave events
if ~isfield(cfg.sw, 'autoFind'), cfg.sw.autoFind = false; end 
if ~isfield(cfg.sw, 'autoFind_fcn'),      cfg.sw.autoFind_fcn = @d2d_detectSlowWaves; end 
% if ~isfield(cfg.sw, 'autoFind_fcn_pars'), cfg.sw.autoFind_fcn_pars = []; end  % default parameters for this are passed later
if ~isfield(cfg.sw, 'ep_lats'),  cfg.sw.peakwin = [-0.05 0.05]; end 

% defaults - ICA
if ~isfield(cfg, 'ica'), cfg.ica = struct; end

%% start timing
tIN_main = tic;

%% load dataset from disk if no variable was passed in from memory  
if isempty(EEG)
    ctemp = cfg;
    if ~any([cfg.fft.run, cfg.sw.run])
        ctemp.loadmode = 'info'; % no need to load actual data unless computing fft across channels or finding slow-waves
    end
    EEG = cleanSleep_loadDataset(ctemp);
end

%% more defaults
if ~isfield(cfg.ica,  'run'), cfg.ica.run  = false;  end

%% get ASR info from the dataset
if EEG.etc.dusk2dawn.asr.raw == false
    fprintf('d2d_validateCleaning: validating ASR results ...\n')    

    % cleaned data  --
    % info & parameters
    valid.raw = false;
    
    % record any parameters that were varied  
    if isfield( EEG.etc.dusk2dawn.asr, 'pars') 
        valid.asr.pars = EEG.etc.dusk2dawn.asr.pars;
    end

    % asr outputs
    valid.asr.propCleaned       = EEG.etc.dusk2dawn.asr.propCleaned;
    valid.asr.propGFPRemoved    = EEG.etc.dusk2dawn.asr.propGFPRemoved;
    valid.asr.propVarRemoved    = EEG.etc.dusk2dawn.asr.propVarRemoved;
    valid.asr.tElapsed          = EEG.etc.dusk2dawn.asr.tElapsed;
    valid.asr.calibDataLen      = EEG.etc.dusk2dawn.asr.calibDataLen;
    
else      
    fprintf('d2d_validateCleaning: using raw data as a reference for ASR validation ...\n')
    % raw data --
    % info & parameters
    valid.raw         = true;
    
    % record any parameters that were varied  
    if isfield( EEG.etc.dusk2dawn.asr, 'pars') 
        valid.asr.pars = EEG.etc.dusk2dawn.asr.pars; % should be nan for raw data
    end

    % asr outputs
    valid.asr.propCleaned     = 0;
    valid.asr.propGFPRemoved  = 0;
    valid.asr.propVarRemoved  = 0;
    valid.asr.tElapsed        = 0;
    valid.asr.calibDataLen    = 0;
    
end

%% compute power spectrum
if cfg.fft.run
fprintf('d2d_validateCleaning: computing frequency decomposition (spectopo) ...\n')
tIN_fft = tic;
fft = cfg.fft; % get sub-structure

% compute spectrum
ents = EEG.event; ents(~strcmp({ents.type},'boundary')) = []; % find boundaries in data e.g. from splitByStage
bounds = [ents.latency];
if ~isempty(bounds)
    [fft.db,fft.freqs,~,~,~] = pop_spectopo(EEG, 1, [], 'EEG', 'plot','off', 'boundaries', bounds );
else
    [fft.db,fft.freqs,~,~,~] = pop_spectopo(EEG, 1, [], 'EEG', 'plot','off');
end
fft.amp = sqrt( 10 .^ (fft.db/10) ); % convert from decibels power to amplitude
fft.amp = fft.amp / length(fft.freqs); % normalise amplitude by N
fft.amp_avgChan   =  mean(fft.amp);

% check if every bin can be computed (i.e. if data are sufficiently sampled)
maxFreq = max(fft.freqs,[],'all');
if maxFreq < max(fft.binFreqs,[],'all')
    warning([mfilename ': data sampling rate is too low for some of the requested frequency bands' ])
end

% bin the FFT means and proportions
nbands = size(fft.binFreqs,1);
fft.binAmp = nan(EEG.nbchan,nbands);
for k = 1 : nbands
    f1 = fft.binFreqs(k,1); f2 = fft.binFreqs(k,2);
    if f2 <= maxFreq
        ind1 = findnearest( fft.freqs , f1 );
        ind2 = findnearest( fft.freqs , f2);
        fft.binAmp(:,k) = mean( fft.amp(:, ind1:ind2) ,2 ); % get amp in bins across channels
    end
end
fft.binAmp_avgChan = mean(fft.binAmp);


%  OLD (custom fft function)
% % compute spectrum
% for c = 1  :EEG.nbchan
%     tempfft = rs_fft(EEG.data(c,:)');
%     if c == 1
%        fft.amp =  nan(size(tempfft,1), size(EEG.chanlocs,2));
%     end
%     fft.amp(:,c) = tempfft;
% %     disp(c)
% end; clear tempfft
% fft.amp = fft.amp';
% fstep = 1/  ((1/EEG.srate)*EEG.pnts); 
% fft.freqs = 0  : fstep : (size(fft.amp,2)-1) * fstep; % construct frequency vector
% 
% % smooth
% smoothLen =  0.5 / fstep; % how many steps needed to smooth at X Hz ?
% if smoothLen > 1
%     fft.amp = movmean(fft.amp, smoothLen, 2, 'Endpoints','fill' ); % mov mean to smooth before downsampling
% end
% 
% % downsample  ? bit weird to recompute the downsample factor each time .. maybe just do by a fixed amount somehow
% for ds = 8:-1:1  
%     if mod( length(fft.freqs) , ds ) == 0 % if divisible by this ds factor
%         break
%     end
% end
% endfreq = findnearest(fft.freqs, fft.binFreqs(end,2)); 
% fft.freqs   =  fft.freqs(1,1:ds:endfreq);
% fft.amp     =  fft.amp(:,1:ds:endfreq); 
% fft.amp_avgChan   =  mean(fft.amp);

% % bin the FFT means and proportions
% for k = 1:size(fft.binFreqs,1)
%     ind1 = findnearest( fft.freqs , fft.binFreqs(k,1) );
%     ind2 = findnearest( fft.freqs , fft.binFreqs(k,2) );
%     fft.binAmp(:,k) = mean( fft.amp(:, ind1:ind2) ,2 ); % get fft in bins across channels
% end
% fft.binAmp_avgChan = mean(fft.binAmp(:,:));
% % remove continuous mean FFT across channels if not requested
% if ~fft.storeWholeFFT
%     fft.freqs        = nan;
%     fft.amp          = nan; 
%     fft.amp_avgChan  = nan;
% end

% outputs
valid.fft = fft;  % store outputs 
tOUT_fft = toc(tIN_fft);
fprintf('d2d_validateCleaning: ... computed frequency decomposition in %.2f mins\n', tOUT_fft/60)
end

%% compute slow-wave (SW) amplitudes
if cfg.sw.run
tIN_sw = tic;
sw = cfg.sw; % get sub-structure

% get slow-waves either from cfg parameters or automatically 
if cfg.sw.autoFind
    fprintf('d2d_validateCleaning: automatically finding slow-wave events using default function (swalldetectnew) ...\n')
    EEG = d2d_detectSlowWaves(EEG, sw.chan, sw.ampThresh);
    sw.codes = {'SW_neg'};
else
    fprintf('d2d_validateCleaning: using user-supplied slow-wave events for validation ...\n')
end

% check for slow-waves provided or detected
if  ~any(   contains( {EEG.event.type},  sw.codes  ) )  
    error('error: d2d_validateCleaning: no slow-wave events provided / none were found by algorithm')
end 

% extract slow-wave epochs to compute sw metrics  
EEG_sw  = pop_select(EEG,'channel' , sw.chan );  
EEG_sw  = pop_epoch(EEG_sw, sw.codes  , sw.peakwin);
sw.waves = squeeze( EEG_sw.data )';   %  figure; plot(EEG_sw.times/1000, mean(sw.waves))
[~,~,~,temp] = ttest( sw.waves );  % figure; plot(EEG_sw.times/1000, sw_tval)
sw.tval = temp.tstat; clear temp

% compute metrics
sw.nwaves           = EEG_sw.trials;
sw.amp_mean         = mean(sw.waves,   [1,2]); 
sw.amp_median       = median(sw.waves, [1,2]);
sw.tval_median      = median(sw.tval);
sw.waves_median     = median(sw.waves, 2); % i.e. median in window around peak for each individual SW

% remove all waveforms to save storage space
sw = rmfield(sw,'waves');

% outputs
valid.sw = sw;  % store outputs
tOUT_sw = toc(tIN_sw);
fprintf('d2d_validateCleaning: ... computed slow-wave metrics in %.2f mins\n\n', tOUT_sw/60)
end

%% compute ICA quality with ICLabel
if cfg.ica.run
fprintf('d2d_validateCleaning: computing ICA quality with ICLabel ...\n')
tIN_ica = tic;
ica = cfg.ica; % get sub-structure

% update ICA results stored in the main part of the EEGlab structure
EEG.icaweights  =   EEG.etc.dusk2dawn.ica.weights;
EEG.icasphere   =   EEG.etc.dusk2dawn.ica.sphere;
if size(EEG.icaweights,1) == size(EEG.icaweights,2) 
    EEG.icawinv     =  inv( EEG.icaweights * EEG.icasphere   );
else
    EEG.icawinv     =  pinv( EEG.icaweights * EEG.icasphere   );
end

% get IC info
ICs = EEG.icaweights * EEG.icasphere * EEG.data;
ica.compvars = sum(EEG.icawinv.^2).*sum((ICs').^2)/((size(EEG.data,1)*size(EEG.data,2))-1); clear ICs % compute component variances
ica.numIC = size(EEG.icaweights,1);
ica.tElapsed  = EEG.etc.dusk2dawn.ica.tElapsed;

% use IClabel to assess ICA quality 
evalc([ 'EEG = iclabel(EEG);' ]);
ica.iclabel.classProb.classLabs   = EEG.etc.ic_classification.ICLabel.classes;
ica.iclabel.classProb.probs       = EEG.etc.ic_classification.ICLabel.classifications;
ica.iclabel.classProb.mean        =   mean( ica.iclabel.classProb.probs );
ica.iclabel.classProb.medi        = median( ica.iclabel.classProb.probs );

%  class probabilities according to my system (4 categories) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ica.iclabel.classProb4.classLabs      = { 'brain', 'bioArtifact', 'recArtifact', 'unknown'   };
ica.iclabel.classProb4.probs(:,1) =     ica.iclabel.classProb.probs(:,1);
ica.iclabel.classProb4.probs(:,2) = sum(ica.iclabel.classProb.probs(:,2:4), 2);
ica.iclabel.classProb4.probs(:,3) = sum(ica.iclabel.classProb.probs(:,5:6), 2);
ica.iclabel.classProb4.probs(:,4) =     ica.iclabel.classProb.probs(:,7);
ica.iclabel.classProb4.mean  =   mean( ica.iclabel.classProb4.probs );
ica.iclabel.classProb4.medi  = median( ica.iclabel.classProb4.probs );

% compute proportion of components reflecting each category (4 categories)
ica.iclabel.classProb4.catNum = zeros( 1 , 4  );  % inefficient but simpler on main script to put in the loop like this (maybe not though)
[maxProb, maxCat] = max( ica.iclabel.classProb4.probs ,[], 2 );
[temp(:,2),temp(:,1)] = hist(maxCat(:),unique(maxCat(:))); %#ok<HIST>
for k = 1:size(temp,1)
    ica.iclabel.classProb4.catNum(temp(k,1)  ) = temp(k,2);
end; clear temp
ica.iclabel.classProb4.catNumProp =  100 * ica.iclabel.classProb4.catNum /  ica.numIC;

% compute summed variance of each category (4 categories)
for k = 1:4
   ica.iclabel.classProb4.catVar(k) =  sum( ica.compvars(maxCat == k) );
end; clear temp
ica.iclabel.classProb4.catVarProp = 100 * ica.iclabel.classProb4.catVar / sum(ica.iclabel.classProb4.catVar);

% compute identifiability difference (diff prob between most likely category and unknown (other) categories)
identiDiffs = max( ica.iclabel.classProb4.probs(:,1:3),[],2)  -  ica.iclabel.classProb4.probs(:,4);
ica.iclabel.identiDiff.allDiffs  = identiDiffs ;
ica.iclabel.identiDiff.meanDiff  = mean(  identiDiffs  );

% compute proportion of components that are most likely identifiable 
knownComps   = sum( maxCat == 1 | maxCat == 2 | maxCat == 3 );
unknownComps = sum( maxCat == 4 );
knownProp   =  100 * knownComps / (knownComps+unknownComps);
ica.iclabel.identiDiff.propNumKnown = knownProp;

% compute proportion of total comp variance accounted for by identifiable components
knownVar    = sum(  ica.compvars(maxCat == 1 | maxCat == 2 | maxCat == 3)   );
unknownVar  = sum(  ica.compvars(maxCat == 4)   );
ica.iclabel.identiDiff.propVarKnown = 100 * knownVar / (knownVar+unknownVar);

% outputs
valid.ica = ica; % store outputs
tOUT_ica = toc(tIN_ica);
fprintf('d2d_validateCleaning: ... computed ICA quality in %.2f mins\n', tOUT_ica/60)
end

%% store outputs in dataset
EEG.etc.dusk2dawn.valid = valid;

%% resave header with validation metrics stored in .etc field
% ctemp = []; ctemp.savePath = cfg.savePath; ctemp.headerOnly = true;
% cleanSleep_saveDataset(EEG,ctemp);

%% timing end
tOUT_main = toc(tIN_main);
fprintf('d2d_validateCleaning: ... file %s finished in %.2f mins\n', EEG.filename, tOUT_main/60)


end
