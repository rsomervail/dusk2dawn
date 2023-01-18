%
%       wrapper function for the vis_artifacts function in the clean_rawdata plugin for EEGLAB
% 
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net
% History:
% 11/07/2022 ver 0.1. Created
% 
%%  

function EEG = d2d_vis_artifacts( EEG, cfg )

%% check whether 2 datasets were selected, if no then just assume 0 varied parameters and plot data before/after ASR
if isempty(cfg)

    % EEG plotted in red
    ctemp = [];
    ctemp.loadRaw = true;
    EEG_red = d2d_applyCleaning(EEG,ctemp);

    % EEG plotted in blue
    ctemp = [];
    ctemp.loadRaw = false;
    EEG_blue = d2d_applyCleaning(EEG,ctemp);

%% else load 2 datasets specified by the cfg selections
else

    % EEG plotted in red
    ctemp = [];
    flds = fields(cfg);
    for k = 1:length(flds) % get parameters needed to load dataset
        ctemp.( strrep(flds{k},'selL','sel') ) = cfg.(flds{k}); % rename field so that lower level function recognises the inputs
    end
    EEG_red = d2d_applyCleaning(EEG,ctemp);
    
    % EEG plotted in blue
    ctemp = [];
    flds = fields(cfg);
    for k = 1:length(flds) % get parameters needed to load dataset
        ctemp.( strrep(flds{k},'selR','sel') ) = cfg.(flds{k}); % rename field so that lower level function recognises the inputs
    end
    ctemp.loadRaw = false;
    EEG_blue = d2d_applyCleaning(EEG,ctemp);

end

%% call vis_artifacts function with these two datasets
vis_artifacts( EEG_blue, EEG_red ); 
