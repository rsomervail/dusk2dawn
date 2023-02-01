%
%  - Saves EEG data processed by Dusk2Dawn, can either save whole dataset or just update the header.
%  
%  
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net 
%%  
function EEG = d2d_saveDataset(EEG,cfg)

% defaults
if ~isfield(cfg,'saveFlag'),   cfg.saveFlag   = true;  end
if ~isfield(cfg,'headerOnly'), cfg.headerOnly = false; end

% process inputs
EEG.setname  = strrep( cfg.saveName, '.set','' ); % ? may be unnecessary as setname should not contain file extension
EEG.filename = [EEG.setname '.set']; 
EEG.filepath = cfg.savePath;

% save dataset
tic
if cfg.saveFlag
    if ~cfg.headerOnly  
        pop_saveset(EEG, 'filepath', EEG.filepath, 'filename', EEG.filename, 'savemode','twofiles' ); 
    elseif cfg.headerOnly
        EEG.data = EEG.datfile; save([ EEG.filepath filesep EEG.filename  ],'EEG') % saving only header because data does not change
    end
    fprintf('... saved file ''%s'' to:\n  --> "%s"\n', EEG.setname, [EEG.filepath filesep EEG.filename]) 
else
    fprintf('... not saving file ''%s'' to disk\n', EEG.setname) 
end
toc

end


