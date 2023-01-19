%
%       apply ASR cleaning using d2d_applyCleaning
% 
%
% 
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net
% History:
% 19/01/2023 ver 1.0.0 Created
% 
%%  
function EEG = d2d_applyCleaning( varargin )

% get inputs
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

% get raw data filename & clean filenames for later
origFile = strrep(EEG.etc.dusk2dawn.cfg.origFile,'.set','');
cleanFiles = EEG.etc.dusk2dawn.cfg.cleanFiles;

% load either raw data or a cleaned file
if cfg.loadRaw

    fprintf('dusk2dawn: applying ASR cleaning; load raw data ...\n')
    file2load = origFile;

else % if loading cleaned data
    
    % find clean file
    fprintf('dusk2dawn: applying ASR cleaning; previously-cleaned dataset will be loaded from disk ...\n')
    if isfield(cfg,'sel_par_1')
        load1 = [ '_p1-' num2str(cfg.sel_par_1,'%02d') ];
    else
        load1 = '';
    end
    if isfield(cfg,'sel_par_2')
        load2 = [ '_p2-' num2str(cfg.sel_par_2,'%02d') ];
    else
        load2 = '';
    end
    if isfield(cfg,'sel_par_3')
        load3 = [ '_p3-' num2str(cfg.sel_par_3,'%02d') ];
    else 
        load3 = '';
    end
    file2load = [ origFile '_clean' load1 load2 load3  ];

end % if loadRaw

% load dataset
fprintf('dusk2dawn: loading dataset : %s.set ...\n', file2load )
EEG = pop_loadset( 'filepath',EEG.filepath,'filename', [file2load '.set'] );

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
