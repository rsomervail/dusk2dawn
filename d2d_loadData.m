%
%       load a version of the current dataset that has been processed with Dusk2Dawn
%           - e.g. apply ASR cleaning with a particular set of parameters
%              or revert back to the raw data
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net 
%%  
function EEG = d2d_loadData( varargin )


%% get inputs
EEG = varargin{1};
if nargin == 1
    cfg = [];
elseif nargin == 2
    cfg = varargin{2};
else
    error(['Error: ' mfilename ': too many input arguments'])
end

% defaults 
if ~isfield(cfg,'loadRaw'), cfg.loadRaw = false; end

% get dataset filenames and paths
savePath = EEG.etc.dusk2dawn.cfg.savePath;
origFile = strrep(EEG.etc.dusk2dawn.cfg.origFile,'.set','');
cleanFiles = EEG.etc.dusk2dawn.cfg.cleanFiles;

% load either raw data or a cleaned file
if cfg.loadRaw

    fprintf('d2d_loadData: loading raw data from disk ...\n')
    file2load = origFile;

else % if loading cleaned data
    
    % find clean file
    fprintf('d2d_loadData: loading a previously-cleaned dataset from disk ...\n')
    file2load = d2d_getCleanedFileName(origFile,cfg);

end % if loadRaw

% load dataset
fprintf('d2d_loadData: dataset = %s.set ...\n', file2load )
EEG = pop_loadset( 'filepath', savePath,'filename', [file2load '.set'] );

% double-check that data has requested parameter values
if ~cfg.loadRaw
    pars = EEG.etc.dusk2dawn.cfg.pars;
    npars = length(pars.values);
    if npars > 0
        checkFlag = true;
        for k = 1:npars
                thisSelection = pars.values{k}(cfg.(['sel_par_' num2str(k)]));
                thisDataset   = EEG.etc.dusk2dawn.valid(1).asr.pars.(['par_' pars.labels{k}]);
                if  ~isequal( thisSelection, thisDataset )
                    checkFlag = false;
                end
        end
    
        % if dataset doesn't match the requested parameter values then give error
        if checkFlag == false
            error( ['Error: ' mfilename ': loaded dataset doesn''t match requested dataset, submit bugreport to: www.github.com/rsomervail/dusk2dawn' ])
            %! if this ever happens then write a subfunction to return all info structs then load the dataset with matching param vals
        end
    end
end


end % END FUNCTION
