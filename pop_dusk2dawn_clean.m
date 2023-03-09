%
%   Usage:
%   >> EEG = pop_dusk2dawn(EEG); % launch a GUI
% 
%   Clean dataset with 1 or more sets of ASR parameters & compute validation metrics
%   (note that output EEG structure is the raw data, and cleaned data must be loaded 
%    with d2d_loadData)
% 
% ## HOW TO USE ##
% 
%   - First you must choose the ASR parameters for the cleaning. The most important parameter here
%     is the cutoff value. This determines the ASR severity, i.e. how much larger artifacts must be 
%     than signal to be removed by ASR, and is expressed in standard deviations (SD). In sleep we 
%     have found values from 30 to 45 SD are appropriate to avoid removing Slow Waves.
%           
%         - You can vary up to 3 parameters in this screen by entering a range of values like this: 
%               [ 30, 35, 40, 45 ]
%               
%         - Advanced settings can also be found by clicking the popup buttons, but for most users 
%           this is unnecessary.
%            
%   - Second, you must choose the ASR implementation. By default the data will be segmented by sleep 
%     stage and ASR performed in a sliding-window. You must specify the events in your dataset which 
%     correspond to sleep stages using the selection box. You can also vary the sliding window length 
%     and overlap, although the default settings are probably fine here.
% 
%   - Third you can choose which additional validation metrics you would like to compute before & after 
%     cleaning. These include:
% 
%         (1) How much variance ASR removed from the data (computed by default)
%                
%         (2) Power by frequency band
%                
%         (3) Slow-wave amplitude and consistency (here you can either select previously-loaded events 
%             in your dataset which correspond to SW events or automatically detect them in each dataset)
%                
%         (4) Quality of Independent Component Analysis (ICA) decomposition using IClabel. This allows you
%             to judge how much ASR has removed artifacts which are difficult to classify using ICA (labelled
%             'unknown' components) and how much of the data is comprised of brain activity.
% 
% 
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net
%%  

function [EEG, cancelFlag ] = pop_dusk2dawn_clean(EEG)
%% HANDLE MULTIPLE FILES
nfiles = length(EEG);
if nfiles > 1
    
    % check that datasets are compatible, they should have the same channels at least (for ICA & SW detection)
    if ~isequal( EEG.chanlocs ) || ~isequal(  EEG.nbchan  )
        errordlg2('! Channels of datasets do not match, this will probably cause errors later if you want to compute slow wave amplitudes or run ICA' ...
            ,'Warning');
    end
    if ~isequal( EEG.srate )
        errordlg2('! Sample-rates of datasets do not match, this will probably cause compatibility issues later' ...
            ,'Warning - Sample-rates of datasets do not match');
    end

    % get default savePath
    if sum( ~cellfun( @isempty, {EEG.filepath} ) ) == nfiles
        if isequal( EEG.filepath )
            savePath_default = EEG(1).filepath;
        else
            savePath_default = cd;
        end
    else
        savePath_default = cd;
    end
    if ~strcmp(savePath_default(end),filesep)
        savePath_default = [savePath_default filesep];
    end
%     saveName_default = 'D2D '; 
    
    % get channel list
    misc_valid.chans = {EEG(1).chanlocs.labels};
    misc_valid.chanlist = cellfun( @(a,b) [ a '_' b ] , strsplit(num2str( [1:EEG(1).nbchan] )), misc_valid.chans , 'UniformOutput', false);
    
    % get number of chans and sample rate for calc of default ASR window
    misc_asr.nbchan = EEG(1).nbchan;
    misc_asr.srate  = EEG(1).srate;

else % if only one dataset

    % get default savePath
    if isempty( EEG.filepath )
        savePath_default = cd;
    else
        savePath_default = [EEG.filepath];
    end
    if ~strcmp(savePath_default(end),filesep)
        savePath_default = [savePath_default filesep];
    end
%     if ~isempty(EEG.filename)
%         saveName_default = strrep( EEG.filename ,'.set','');
%     else
%         saveName_default = strrep( EEG.setname ,'.set','');
%     end
    
    % get channel list
    misc_valid.chans = {EEG.chanlocs.labels};
    misc_valid.chanlist = cellfun( @(a,b) [ a '_' b ] , strsplit(num2str( [1:EEG.nbchan] )), misc_valid.chans , 'UniformOutput', false);
    
    % get number of chans and sample rate for calc of default ASR window
    misc_asr.nbchan = EEG.nbchan;
    misc_asr.srate  = EEG.srate;

end


%% set defaults
% handle core function defaults
sw_autoFind_default  = false;
splitByStage_default = false;

% handle advanced settings defaults
cfg_valid = advanced_popup_valid([],false, misc_valid); % this returns defaults for the advanced validation options popup
cfg_asr   = advanced_popup_asr([],false, misc_asr); % this returns defaults for the advanced asr options popup

%% get relevant events (stages & slow-waves) & set default selections
if sum(  ~arrayfun(  @(x) isempty(x.event), EEG  ) ) == nfiles % if any EEG input has no events, then there will ofc be no common events

    % get all shared event codes across all datasets selected
    allEnts = arrayfun(  @(x) unique({x.event.type}), EEG  , 'UniformOutput', false);
    allEnts = [allEnts{:}];
    if nfiles > 1
        uni = unique(allEnts);
        indy = cellfun( @(x) sum(strcmp(x,allEnts)) ,  uni  );
        ents = uni( indy==nfiles ); % find event codes which are present in all datasets
    else
        ents = allEnts;
    end
    ents = [ ' ' ents ]; % append with nothing so that it is not obligatory to make a selction

    % find any that look like slow waves and set as defaults
    ents_sw_default = find(  contains(ents,'sw','IgnoreCase',true)); %#ok<*EFIND> 
    ents_st_default = find( ~contains(ents,'sw','IgnoreCase',true) & ~strcmp(ents,' ') );
    if isempty(ents_sw_default)
        sw_autoFind_default  = true; % if no clear SW events exist then default is to automatically detect them
    end
    if ~isempty(ents_st_default) 
        splitByStage_default = true; % if possible stage events exist, default is to split by them
    end

else % if no events found at all 
    ents_sw_default = [];
    ents_st_default = [];
    sw_autoFind_default  = true; % if no events found at all, default is to auto-find SWs
    splitByStage_default = false; % if no events then can't split by stage
end

% % set to 1, i.e. blank value
if isempty(ents_sw_default), ents_sw_default = 1; end
if isempty(ents_st_default), ents_st_default = 1; end

%% DEFINITIONS

% geo (horizontal)
% indent1 = 0.05; % default
indent1 = 0.1;
button_len = [1 3]; % horizontal length of buttons is ratio of first to second numbers

% UI
space = { {} }; % empty line

% vertical
vtitle1 = 1.2;
v0 = 0.2; % empty line height - very small
v1 = 0.5; % empty line height - small
v2 = 0.8; % empty line height - big
veditbox1 = 1; % edit box height
veditbox2 = 1.2; % edit box height
vbutton = 1; % button height

%% DEFINE GUI ELEMENTS - LINE __________
linelen   = 140;
geo_line  = { 1 };
vert_line = [ 1 ];
ui_line   = { { 'Style', 'text', 'string', repmat('_',1,linelen), 'fontweight','bold' ,'fontsize',13} };

%% DEFINE GUI ELEMENTS - GENERAL OPTIONS


cmd = [ ' '' d2d_getSettingsFromSet; '' ' ]; %#ok<NBRAK2> 
callback_getset = [' evalin(''caller'', ' cmd ' )   ']; % callback for button which gets settings from previously-cleaned file

geo_general = { ...
    1   ... % WELCOME TEXT 1
    1   ... % WELCOME TEXT 2
    [3 0.5  1.1] ... % savePath editbox
...    [1 0.4  2.3]   ... % saveName editbox 
};
vert_general = [ ...
    v2   ... % WELCOME TEXT 1
    v2   ... % WELCOME TEXT 2
    veditbox2 ... % savePath editbox + title
...     veditbox2   ... % saveName editbox + title
];
ui_general = { ...
    { 'Style', 'text', 'string', 'Welcome to Dusk2Dawn!   -  you can choose below how to clean your data and which metrics to compute for later validation of the cleaning results', 'fontweight', 'bold','fontsize',13, 'horizontalalignment', 'Center' } ...
    { 'Style', 'text', 'string', [repmat(' ',1,12) '>> You can also enter multiple values for up to 3 options and test the effects of each on your data, e.g. asr_cutoff = [30,35,40,45]'], 'fontangle', 'italic', 'horizontalalignment', 'center' } ...
    { 'Style', 'edit', 'string', savePath_default 'tag' 'savePath' } ...   
        {} ... ... { 'Style', 'text', 'string', '<--- folder' 'fontweight', 'bold' 'fontangle','italic' } ...
        { 'Style','pushbutton', 'string','Load options from set file', 'callback',callback_getset} ...
... %     { 'Style', 'edit', 'string', saveName_default 'tag' 'saveName' } ...
... %         {} ... { 'Style', 'text', 'string', '<--- save as' 'fontweight', 'bold' 'fontangle','italic'} ...
... %         {} ...
};

%% DEFINE GUI ELEMENTS - ASR OPTIONS - ASR PARAMETERS
geo_asr_params = { ...
    1 ... % ASR OPTIONS
    1 ... % 
    [1 3] ... % asr cutoff editbox + text 
    [1] ... % help text line 1
};
vert_asr_params = [ ...
    vtitle1 ... % ASR OPTIONS
    v0 ... %
    veditbox1 ... % asr cutoff editbox
    1 ... % help text
];
ui_asr_params = { ...
    { 'Style', 'text', 'string', ['ASR OPTIONS  ' repmat('_',1,linelen+17)], 'fontweight', 'bold', 'fontsize',13} ...
    {} ...
    { 'Style', 'edit', 'string', '[  35  ]' 'tag' 'asr_cutoff'  } ... 
        { 'Style', 'text', 'string', '<--- enter ASR cutoff(s)' 'fontweight', 'bold'} ...
    { 'Style', 'text', 'string', '      (determines ASR severity, lower cutoff = more severe filtering - we recommend trying several values, e.g. [30,35,40,45])', 'fontangle','italic'} ...  
};

%% DEFINE GUI ELEMENTS - D2D OPTIONS
geo_d2d = { ... 
    1 ... % split by stage checkbox
    [ indent1  .2 1 ]   ...  % stage codes listbox + instruction
    1 ...
    [ indent1  .2 1 ] ... % stage window editbox + instruction
    1 ... % run ASR in a sliding window?
    [ indent1 .2 1] ...  % sliding window length editbox + instruction
};
vert_d2d = [ ...
    1 ... % split by sleep stage before ASR cleaning?
    length(ents) ... % stage codes listbox + instruction
    v0 ...
    veditbox1 ... % stage window editbox + instruction
    1 ... % run ASR in a sliding window?
    veditbox1 ... % sliding window length editbox + instruction
];
ui_d2d = { ...
    { 'Style', 'checkbox', 'string' ' (1) Split data by sleep stage before ASR cleaning?' 'fontweight', 'bold' 'value' splitByStage_default 'tag' 'splitByStage' } ...
    {} ...
        { 'Style','listbox','string',ents,'value',ents_st_default,'Max',length(ents),'tag','stageCodes'} ... 
        { 'Style', 'text', 'string', '<--- select event codes which mark sleep stages' 'fontweight', 'bold' } ...
    {} ...
    {} ...
        { 'Style', 'edit', 'string', '[ 0, 30 ]' ,  'tag' 'stageWin'  } ...
        { 'Style', 'text', 'string', '<--- choose window around sleep stage events (s)' 'fontweight', 'bold'  } ...
    { 'Style', 'checkbox', 'string' ' (2) Run ASR in a sliding window?' 'fontweight', 'bold' 'value' ~splitByStage_default 'tag' 'splitBySlidingWindow' } ...    
    {} ...    
        { 'Style', 'edit', 'string', '[   8   ]' 'tag' 'chunk_len'  } ...    
        { 'Style', 'text', 'string', '<--- choose length of sliding window (mins)','fontweight', 'bold'} ...
};

%% DEFINE GUI ELEMENTS - ADVANCED ASR/D2D SETTINGS

cmd = [ ' '' cfg_asr = advanced_popup_asr(cfg_asr,true, misc_asr); '' ' ];
callback = [' evalin(''caller'', ' cmd ' )   '];
geo_adv_asr =  { ...
    button_len ...
};
vert_adv_asr = vbutton * ones(size(geo_adv_asr));  
ui_adv_asr = { ...
    { 'Style','pushbutton', 'string','Advanced ASR Options', 'callback',callback } ...
        { 'Style', 'text', 'string', '(generally you can leave these options unchanged)','fontangle','italic'} ...
};



%% DEFINE GUI ELEMENTS - VALIDATION OPTIONS

geo_valid = { ...
        1 ... % POST-CLEANING VALIDATION 
        1 ... % validation - FFT
        1 ... % validation - SW  
        [ indent1  .2 1 ] ...  % SW event selection box + title
        [ indent1  .2  ] ... % SW OR detect automatically checkbox
        1 ... % validation - ICA 
        [1 1 1] ... % ** requires channel locations **
};
vert_valid = [ ...
        vtitle1 ... % POST-CLEANING VALIDATION 
        1 ... % validation - FFT
        1 ... % validation - SW  
        length(ents) ...  % SW event selection box + % SW OR detect automatically checkbox
        v2 ... % SW OR detect automatically checkbox
        1 ... % validation - ICA 
        1 ... % ** requires channel locations **
];
ui_valid = { ...
        { 'Style', 'text', 'string', ['POST-CLEANING VALIDATION  ' repmat('_',1,linelen)],'fontsize',13, 'fontweight', 'bold' } ...
        { 'Style', 'checkbox', 'string' ' (1) Compute frequency spectrum?' 'value' 0 'tag' 'fft_run', 'fontweight', 'bold' } ...
        { 'Style', 'checkbox', 'string' ' (2) Measure slow-wave amplitude?' 'value' 0 'tag' 'sw_run', 'fontweight', 'bold' } ...
        {} ...
            { 'Style','listbox','string',ents,'value',ents_sw_default,'Max',length(ents),'tag','sw_codes'} ...  
            { 'Style', 'text', 'string', '<--- select event codes which mark slow waves','fontweight', 'bold'} ...    
        {} ...
            { 'Style', 'checkbox', 'string' '<--- OR detect slow-waves automatically instead','fontweight', 'bold', 'value' sw_autoFind_default 'tag' 'sw_autoFind' } ...     
        { 'Style', 'checkbox', 'string' ' (3) Assess ICA quality before/after ASR using IClabel?' 'value' 0 'tag' 'ica_run', 'fontweight', 'bold' } ... 
        { 'Style', 'text', 'string', '         ** requires channel locations **','fontangle','italic' } ...
            {} ...
            {} ...
    };

%% DEFINE GUI ELEMENTS - ADVANCED VALIDATION SETTINGS
cmd = [ ' '' cfg_valid = advanced_popup_valid(cfg_valid,true,misc_valid); '' ' ];
callback = [' evalin(''caller'', ' cmd ' )   '];
geo_adv_valid =  { ...
    button_len ...
};
vert_adv_valid = vbutton * ones(size(geo_adv_valid));  
ui_adv_valid = { ...
    { 'Style','pushbutton', 'string','Advanced Validation Options', 'callback',callback } ...
        {} ...
};

%% BUILD GUI 
geometry = [ ...
    geo_general ...
    1 ...
    geo_asr_params ...
    geo_d2d ...
    1 ...
    geo_adv_asr ...
    1 ...
    geo_valid ...
    1 ...
    geo_adv_valid ...
];
vert = [ ...
    vert_general    ...
    v0 ...
    vert_asr_params ...
    vert_d2d        ...
    v0               ...
    vert_adv_asr    ...
    v0 ... 
    vert_valid      ...
    v0               ...
    vert_adv_valid  ...
];
uilist = [ ...
    ui_general ...
    space ...
    ui_asr_params ...
    ui_d2d ...
    space ...
    ui_adv_asr ...
    space ...
    ui_valid ...
    space ...
    ui_adv_valid ...
];


%% RUN GUI AND GET PARAMETERS

[ ~,~,strhalt, cfg ] = inputgui( 'geometry',geometry,'uilist', uilist, 'geomvert',vert, ...
   'helpcom','pophelp(''pop_dusk2dawn'');', ...
   'title','Clean whole-night dataset using ASR -- pop_dusk2dawn_clean()' ...
   ,  'skipline', 'off' );

% return if cancelled
if ~strcmp(strhalt,'retuninginputui') % if cancelled
    cancelFlag = true;
    return
else 
    cancelFlag = false;
end

%% UPDATE ADVANCED SETTINGS
% need to copy advanced settings to cfg here using the subfunction
cfg = updateCFG(cfg,cfg_valid);
cfg = updateCFG(cfg,cfg_asr);

%% make save directory if doesn't exist yet
if ~exist(cfg.savePath,'dir')
    mkdir(cfg.savePath);
end

%% check for matching files in this directory, indicating there is a run of D2D already present
files = dir(cfg.savePath);
files = files(~[files.isdir]);
files = {files.name};           
files(~endsWith(files,'.set')) = []; % remove all but .set files
files2clean = {EEG.filename};
files2clean = strrep(files2clean,'.set',''); % remove extensions
if any(  ...
        contains(files,'_clean') ...
        &  cellfun(  @(x) any(startsWith(x, files2clean )) , files )  ...
      )   
    errordlg2('! looks like there are already D2D-cleaned files in this folder - this may cause compatibility issues' ...
        ,'Error - The chosen save folder may contain a previous D2D runthrough for these datasets');
%     error 'looks like there are already D2D-cleaned files in this folder - delete them or choose another folder';
end

%% store & reformat parameters
% ASR settings
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

% validation settings - fft
cfg.fft.run = cfg.fft_run; cfg = rmfield(cfg,'fft_run');
evalc([ 'cfg.fft.binFreqs = ' cfg.fft_binFreqs ]); cfg = rmfield(cfg,'fft_binFreqs');
evalc([ 'cfg.fft.binFreqsLabels = ' cfg.fft_binFreqsLabels ]); cfg = rmfield(cfg,'fft_binFreqsLabels');

% validation settings - sw
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

% validation settings - ica
cfg.ica.run = cfg.ica_run; cfg = rmfield(cfg,'ica_run');
evalc([ 'cfg.ica.numIC = ' cfg.ica_numIC ]); cfg = rmfield(cfg,'ica_numIC');

%% % check options were appropriate for these data
% if cfg.ica.run && sum(isempty(),'all') % ? not sure which chanlocs are necessary for IClabel
%     warning 'pop_dusk2dawn: '          % ? but EEGLAB has a binary display for whether dataset has locs
% end

% if trying to stage split but no events selected for sleep stages then throw error
if cfg.splitByStage && isempty(cfg.stageCodes)
    errordlg2('Please try again and select the events marking sleep stages, or do not split by sleep stage' ...
        ,'Error: no sleep stage events were selected');
    return
end
% if trying to validate slow waves but no events selected for SWs then throw error
if cfg.sw.run && isempty(cfg.sw.codes) && ~cfg.sw.autoFind
    errordlg2('Please try again and select the events marking slow waves, find them automatically, or do not validate slow waves' ...
        ,'Error: no slow wave events were selected');
    return
end

% if selecting SW events manually, then check SW events and sleep stage events don't overlap (if using both stage split & SW user-provided events)
if ~cfg.sw.autoFind && cfg.sw.run && cfg.splitByStage 
    if any(contains( cfg.stageCodes, cfg.sw.codes )) || any(contains( cfg.sw.codes, cfg.stageCodes ))
        errordlg2('Please try again and make sure the selected slow wave and sleep stage events are different ' ...
            ,'Error: you selected some events as both sleep stages and slow waves');
        return
    end
end

%% call dusk2dawn cleaning function 
EEG = dusk2dawn_clean(EEG, cfg); % <-------------------------------------------------------

%% SUBFUNCTION - update cfg with advanced options

% update config cfg
function cfg = updateCFG(cfg,add)
    flds = fields(add);
    for f = 1:length(flds)
        cfg.(flds{f}) = add.(flds{f});
    end
end


%% END MAIN FUNCTION
end
