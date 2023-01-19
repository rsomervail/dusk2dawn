%
%
%       - takes continuous whole-night sleep EEG data in EEGlab format
%       - checks whether the data are adequate for ASR pipeline:
%           - checks whether data are continuous
%           - checks whether data have enough channels for ASR
%           - checks whether data are average referenced
%           - checks for reference channels still in data
%           - checks for interpolated channels  
% 
% 
%
% 
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net
% History:
% 19/01/2023 ver 1.0.0 Created
% 
%%  
function EEG = d2d_checkData(EEG)

tIN = tic;
fprintf('checking whether data format is adequate for ASR cleaning ...\n')

% verify EEGlab structure
EEG = eeg_checkset( EEG );

% check whether data are continuous
if EEG.trials > 1
   error 'Error: data are epoched -> use only continuous data for this pipeline' 
end

% check whether data have a sufficient number of channels for ASR to run effectively
if      EEG.nbchan < 8
    warning('number of channels is low (%s) -> ASR performance may suffer -> make sure to check the effect of cleaning',num2str(EEG.nbchan))
end

% check if data are average referenced (using EEGlab header information)
% tic
refs = strsplit( strjoin(unique({EEG.chanlocs.ref})) );
refs = refs(cellfun(@(x) ~isempty(x), refs));
if any(contains(refs,'average'))
   error 'Error: data are average-referenced; this is not compatible with the ASR algorithm -> avoid average-referencing before using this pipeline'
end

% check if data are average referenced (by analysing a subset of samples)
samps2check = randperm(EEG.pnts, min(EEG.pnts,1000) ); % get a random sample of the data
if std( mean(EEG.data(:,samps2check)) ) < 0.001
    warning(['mean across channels is very close to zero; data may have been average-referenced prior to import' ...
        '-> if this is not the case then disregard this message, ' ...
        'otherwise avoid average-referencing before using this pipeline '])
end
% toc

% Check for reference channels still in data  
% tic
try
%     tic
    refindy_all = false(1,EEG.nbchan);
    for r = 1:length(refs)
        refindy = strcmp({EEG.chanlocs.labels}, refs{r});
        if any(refindy)
            warning('data contain reference channel: %s',refs{r})
        end
        refindy_all = refindy_all | refindy;
    end
%     toc
    
    % remove any reference channels found
    if any(refindy_all)
        warning('-> removing reference channels, for speed remove these reference channels before using the pipeline')
        EEG = pop_select(EEG,'channel', find(~refindy_all) );
    end
    
catch
    fprint( 'Reference channels are not labelled or channel locations are not present in data\n')
end
% toc

% check for interpolated channels
if contains( EEG.history , 'pop_interp' )
    error 'Error: data contains interpolated channels -> avoid interpolating channels until after ASR'
end

toc(tIN)
fprintf('... data format is adequate\n')


end
