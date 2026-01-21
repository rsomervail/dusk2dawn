% 
%  - split data by sleep stage 
%       - segments of each stage will be concatenated as continuous data, rather than epoched
% 
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net 
%%  
function [EEG_out_all, EEG_dummy ] = d2d_splitByStage(EEG, cfg)

fprintf('splitting data by sleep stage ...\n')
tIN = tic;

% first extract dummy version of the data for later use when recombining the different stages
EEG_dummy = EEG; EEG_dummy.data = nan(size(EEG.data));

% get sleep stages 
if isempty(cfg.stageCodes)
    error('Error: no sleep-stage event codes provide -> please provide these events using the cfg.stageCodes parameter')
end   
fprintf('using the following event codes to separate sleep-stages:\n')
stages = cfg.stageCodes;
disp(stages)

% default parameters
if ~isfield(cfg,'saveStagesToDisk'),   cfg.saveStagesToDisk = false; end
% if ~isfield(cfg,'loadInChunks'),   cfg.loadInChunks = false; end
if ~isfield(cfg,'keepStagesInMemory'), cfg.keepStagesInMemory = true; end
if ~isfield(cfg,'savePaths'), cfg.savePaths = repmat({cd},length(stages),1); end
if ~isfield(cfg,'useEventDurations'), cfg.useEventDurations = false; end
if ~isfield(cfg,'cleanUnscored'), cfg.cleanUnscored = true; end

    % segment stages
    nstages = length(cfg.stageCodes);
    indy2take_all = false(1,EEG.pnts);
    st = 1;
    while st <= nstages
        
        fprintf('processing sleep stage: %s\n',cfg.stageCodes{st})
        
        % find timepoints corresponding to this stage
        if ~strcmp( cfg.stageCodes{st}, 'UNSCORED' ) % check whether this stage is present in events or is the remaining "unstaged" data segments
            ents2take = EEG.event(strcmp( {EEG.event.type}, cfg.stageCodes{st})); % get all events corresponding to this stage
            indy2take = false(1,EEG.pnts); % get all timepoints corresponding to this stage
            if cfg.useEventDurations % if there are durations then just use these
                for e = 1:length(ents2take)
                    indy2take( ents2take(e).latency : ents2take(e).latency+ents2take(e).duration-1 ) = true;
                end
            else % otherwise use intervals parameter to epoch around each sleep stage event
                stageInts = cfg.stageWin*EEG.srate;
                if size(stageInts,1) == 1
                    stageInts = ones(length(ents2take),2).*stageInts;
                end
                for e = 1:length(ents2take)
                    indy2take( ents2take(e).latency+stageInts(e,1) : ents2take(e).latency+stageInts(e,2)-1 ) = true;
                end
            end
        else % if this stage corresponds to the remaining "unscored" segments 
            indy2take = ~indy2take_all;            
        end 

        % extract data from this sleep stage
        evalc( ' EEG_out = eeg_eegrej(EEG, logical2indices(~indy2take) ); '); % negate because function rejects timewindows
        
        % store info about this stage
        EEG_out.etc.dusk2dawn.stageSplit.indy_wholeData = indy2take; % indy of each sample for this stage in whole dataset
            ents = EEG_out.event; ents(~strcmp({ents.type},'boundary')) = [];
        EEG_out.etc.dusk2dawn.stageSplit.bounds = [ents.latency]; % store discontinuities in resulting data after split, for alter recombination
        EEG_out.etc.dusk2dawn.stageSplit.thisStage  = st;
        EEG_out.setname = [EEG_out.setname '_' cfg.stageCodes{st} ];
        
        
        % optionally keep dataset in memory
        if cfg.keepStagesInMemory
            EEG_out_all(st) = EEG_out;
            fprintf('... stored stage ''%s'' in memory\n', cfg.stageCodes{st})
        end
        
        % check for remaining unstaged segments of data
        if cfg.cleanUnscored
            indy2take_all = indy2take_all | indy2take;
            if st == nstages
                if any(~indy2take_all) % check for any remaining data segments 
                    fprintf('identified data segments which have not been scored (or have not been specified by user ...)\n')
                    fprintf('... cleaning these remaining "unscored" data segments ... \n')
%                     fprintf('... cleaning these remaining "unscored" data segments ... (if you don''t want to clean these segments, uncheck "also clean data segments which have not been scored?")\n')
                    % create new stage and continue loop
                    nstages = nstages + 1; 
                    cfg.stageCodes{end+1} = 'UNSCORED';
                end
            end
        end
    
        st = st + 1;
    end % end stage loop 
    for st = 1:nstages
        EEG_out_all(st).etc.dusk2dawn.stageSplit.stageCodes = cfg.stageCodes;
    end
    
% end % end if load chunks 
toc(tIN)
fprintf('... finished splitting data by sleep-stage\n\n')

end % end function
