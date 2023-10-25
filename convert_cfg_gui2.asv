%
%
%   Takes the configuration structure produced by the GUI pop_dusk2dawn_clean GUI (cfg) and
%   reformats the values of each field to make them appropriate for dusk2dawn_clean (cfgout).
%   This will generally not need to be used by any user.
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net
%%
function cfgout = convert_cfg_gui2(cfg,ents)

%% ASR settings
if isempty(ents) 
    cfg.stageCodes = []; 
else
    cfg.stageCodes = ents(cfg.stageCodes);
    cfg.stageCodes(strcmp(cfg.stageCodes,' ')) = []; % remove blank option
    if isempty(cfg.stageCodes), cfg.stageCodes = []; end % make sure it is empty double if empty now
end
evalc([ 'cfg.stageWin = ' cfg.stageWin ]);
if ~cfg.splitBySlidingWindow
    cfg.chunk_len     = nan;
    cfg.chunk_overlap = nan;
else
    evalc([ 'cfg.chunk_len = ' cfg.chunk_len ]);
    evalc([ 'cfg.chunk_overlap = ' cfg.chunk_overlap ]);
end
evalc([ 'cfg.ref_maxbadchannels = ' cfg.ref_maxbadchannels ]);
cfg.ref_maxbadchannels  = cfg.ref_maxbadchannels / 100; % convert MaxBadChannels to proportions for clean_windows function
evalc([ 'cfg.ref_tolerances = ' cfg.ref_tolerances ]);
evalc([ 'cfg.ref_wndlen = ' cfg.ref_wndlen ]);
evalc([ 'cfg.asr_cutoff = ' cfg.asr_cutoff ]);   %#ok<*EVLEQ> 
evalc([ 'cfg.asr_windowlength = ' cfg.asr_windowlength ]);   %#ok<*EVLEQ> 
evalc([ 'cfg.asr_maxdims = ' cfg.asr_maxdims ]);   %#ok<*EVLEQ> 
evalc([ 'cfg.asr_UseRiemannian = ' cfg.asr_UseRiemannian ]);   %#ok<*EVLEQ> 
evalc([ 'cfg.asr_MaxMem = ' cfg.asr_MaxMem ]);   %#ok<*EVLEQ> 
evalc([ 'cfg.maxreftime = ' cfg.maxreftime ]);   %#ok<*EVLEQ> 

%% validation settings - fft
cfg.fft.run = cfg.fft_run; cfg = rmfield(cfg,'fft_run');
evalc([ 'cfg.fft.binFreqs = ' cfg.fft_binFreqs ]); cfg = rmfield(cfg,'fft_binFreqs');
evalc([ 'cfg.fft.binFreqsLabels = ' cfg.fft_binFreqsLabels ]); cfg = rmfield(cfg,'fft_binFreqsLabels');

%% validation settings - sw
cfg.sw.run = cfg.sw_run; cfg = rmfield(cfg,'sw_run');
if isempty(ents) 
    cfg.sw.codes = []; 
else
    cfg.sw.codes = ents(cfg.sw_codes);
    cfg.sw.codes(strcmp(cfg.sw.codes,' ')) = []; % remove blank option if selected
    if isempty(cfg.sw.codes), cfg.sw.codes = []; end % make sure it is empty double if empty now
end
cfg = rmfield(cfg,'sw_codes');
evalc([ 'cfg.sw.peakwin = ' cfg.sw_peakwin ]); cfg = rmfield(cfg,'sw_peakwin');
cfg.sw.chan = cfg.sw_chan; cfg = rmfield(cfg,'sw_chan');
cfg.sw.autoFind = cfg.sw_autoFind; cfg = rmfield(cfg,'sw_autoFind');
evalc([ 'cfg.sw.ampThresh = ' cfg.sw_ampThresh ]); cfg = rmfield(cfg,'sw_ampThresh');
cfg.sw.peakwin = cfg.sw.peakwin / 1000; % convert from ms to s

%% validation settings - ica
cfg.ica.run = cfg.ica_run; cfg = rmfield(cfg,'ica_run');
evalc([ 'cfg.ica.numIC = ' cfg.ica_numIC ]); cfg = rmfield(cfg,'ica_numIC');

%% output
cfgout = cfg;

end