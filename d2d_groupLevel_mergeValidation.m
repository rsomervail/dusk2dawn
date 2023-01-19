% 
% merge validation metrics at group-level
% 
%
% 
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net
% History:
% 19/01/2023 ver 1.0.0 Created
% 
%%  
function cfg = d2d_groupLevel_mergeValidation(cfg)

%% defaults
if ~isfield(cfg,'savePath'),    cfg.savePath  = []; end
if ~isfield(cfg,'loadPaths'),   cfg.loadPaths = []; end

%% if no load or save paths provided, get these with GUI

% get load paths
if  isempty(cfg.loadPaths) % if no data or load paths are passed in at all, use uigetfile to provide paths
   paths = uigetdir2(cd, 'select folders containing the files whose metrics you want to merge');
   cfg.loadPaths = cellfun( @(x) [x '\EEG_orig.set'], paths, 'UniformOutput', false);
end

if isempty(cfg.savePath)
    [file, path] = uiputfile('.mat', 'save group-level validation metrics');
    cfg.savePath = [path '\' file];
end

%% loop through files and load cross-parameter validation metrics
nfiles = length(cfg.loadPaths); 
for f = 1:nfiles
    
    ctemp = [];
    ctemp.loadPath = cfg.loadPaths{f};
    ctemp.loadMode = 'info';
    EEG = cleanSleep_loadDataset(ctemp);
    v(f) = EEG.etc.cleanSleep.valid2;
    
    % check whether the varied parameters have been provided, in order to overwrite the existing ones
    if isfield(cfg,'pars')
        v(f).pars = cfg.pars;
    end
    
end 

% check that parameters are identical across files and store if identical
if ~isequaln(v.pars)
    error 'cleanSleep_groupLevel_merge: Error: parameters are not identical across files that are being merged'
else
    g.pars = v(1).pars; % store in group metrics if identical
end

%% merge across files
flds = fields(v); flds = flds(~strcmp(flds,'pars')); % get all fields except pars

% loop through subfields containing metrics
for f = 1:length(flds)
    
    % get subfield
    tempIN = [v.(flds{f})]; 
    clear tempOUT
    
    % loop through metrics and concatenate across files
    mets = fields(tempIN); mets = mets(startsWith(mets,'m_')); % get list of metrics
    for m = 1:length(mets)
        tempmet = cat( ndims([tempIN.(mets{m})])+1,  tempIN.(mets{m}) ); % get this metric, concatenated across final dimension + 1
        tempmet = squeeze(tempmet); % squeeze to remove singleton dimensions if there are any single-value ASR parameters
        tempOUT.([ 'g_'  mets{m}(3:end)]) = permute(tempmet,[ndims(tempmet), 1:ndims(tempmet)-1]); % permute so that concat dimension comes first
    end
    
    % get other info related to these datasets and store
    otherInfo = fields(tempIN); otherInfo = otherInfo(~startsWith(otherInfo,'m_')); 
    if ~isempty(otherInfo)
        for o = 1:length(otherInfo)
            if ~startsWith(otherInfo{o},'waves')
                tempOUT.(otherInfo{o}) = tempIN.(otherInfo{o});
            end
        end
    end
    
    % store subfield
    g.(flds{f}) = tempOUT; 
       
end

%% outputs
% store dataset names in group-level validation structure
files = cfg.loadPaths;
for f = 1:length(files)
    temp = strsplit(files{f},{'\','/'});
    g.datasets{f} = temp{end-1};
end

% also store the name used to save this set of group-level metrics
temp = strsplit(cfg.savePath,{'\','/'}); 
g.name = temp{end}(1:end-4);

% store group-level validation structure in cfg
cfg.validG = g;

% save validation file
validG = g;
save([cfg.savePath ],'validG')
fprintf('cleanSleep_groupLevel_merge: metrics merged across datasets\n')
fprintf('saved to:\n')
disp([cfg.savePath])

end % end function
