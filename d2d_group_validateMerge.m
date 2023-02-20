% 
%   merge results of validation from several different datasets at group-level
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net
function EEG = d2d_group_validateMerge(EEG)

%% get validation structure for each dataset
nfiles = length(EEG);
if nfiles == 1
    errordlg2('Please try again first selecting multiple datasets: Datasets -> Select Multiple Datasets ' ...
        ,'Error: you have only selected one dataset');
    return 
end

% check datasets are compatible with each other
d2d_checkSameRun(EEG);

% merge validation structures across datasets
for f = 1:nfiles
    v(f,:) = EEG(f).etc.dusk2dawn.valid_merged;
end

fprintf('%s: merging validation results of the following datasets: \n',mfilename)
disp({EEG.filename}')


%% merge validation structure across datasets

% get all fields which must be merged
flds = fields(v(1)); 
flds = flds(~strcmp(flds,'pars')); 
flds = flds(~strcmp(flds,'stage')); 

% loop through stages
if EEG(1).etc.dusk2dawn.cfg.splitByStage
    nstages = length(EEG(1).etc.dusk2dawn.cfg.stageCodes);
else
    nstages = 1;
end
for st = 1:nstages
 
    % loop through subfields containing metrics
    for f = 1:length(flds)
        
        % get subfield
        tempIN = [v(:,st).(flds{f})]; 
        clear tempOUT
        
        % loop through metrics and concatenate across files
        mets = fields(tempIN); mets = mets(startsWith(mets,'m_')); % get list of metrics
        for m = 1:length(mets)
            tempmet = cat( ndims([tempIN.(mets{m})])+1,  tempIN.(mets{m}) ); % get this metric, concatenated across final dimension + 1
            tempmet = squeeze(tempmet); % squeeze to remove singleton dimensions if there are any single-value ASR parameters
            tempOUT.([ mets{m} ]) = permute(tempmet,[ndims(tempmet), 1:ndims(tempmet)-1]); % permute so that concat dimension comes first
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
        g(st).(flds{f}) = tempOUT; 
           
    end % field

end % stage

%% output groupLevel validation
datasets = {EEG.filename};
for f = 1:nfiles
    EEG(f).etc.dusk2dawn.group.validG = g;
    if ~isfield(EEG(f).etc.dusk2dawn.group,'datasets') % if datasets haven't been merged yet then add filenames about all datasets in group
        EEG(f).etc.dusk2dawn.group.datasets = datasets;
    end
end

fprintf('%s: ... finished merging\n',mfilename)


end % function