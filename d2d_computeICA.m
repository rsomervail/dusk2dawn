%
%
%   - wrapper function for running ICA in the Dusk2Dawn toolbox
%
%
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net 
%%  
function [EEG_out, cfg] = d2d_computeICA(EEG, cfg)

% defaults
if ~isfield(cfg,'ica'),  cfg.ica  = struct; end
% if ~isfield(cfg,'pica'), cfg.pica = struct; end
if ~isfield(cfg.ica, 'run'),       cfg.ica.run        = true;   end
if ~isfield(cfg.ica, 'numIC'),     cfg.ica.numIC      = [];     end
% if ~isfield(cfg.pica,'run'),       cfg.pica.run       = false;  end
% if ~isfield(cfg.pica,'algorithm'), cfg.pica.algorithm = 'lap';  end

%% load dataset from disk if no variable was passed in from memory  
if isempty(EEG)
    EEG = cleanSleep_loadDataset(cfg);
end

%% clean any existing weights in the dataset
EEG.icasphere  = [];
EEG.icaweights = [];
EEG.icawinv    = [];
EEG.icaact     = [];
EEG.etc.oldicachansind = {};
EEG.etc.oldicasphere   = {};
EEG.etc.oldicaweights  = {};

%% compute ICA
% if cfg.ica.run
    fprintf('running ICA ...\n')
    tin_ica = tic;

    % if no number of comps specified just do square ICA
    if isempty(cfg.ica.numIC)
        cfg.ica.numIC = EEG.nbchan;
    end

    % run ICA
    if cfg.ica.numIC == EEG.nbchan
        [EEG] = pop_runica(EEG, 'icatype', 'runica', 'interrupt','off' );  %  ,'stop', 1E-7  );
    else
        [EEG] = pop_runica(EEG, 'icatype', 'runica', 'pca', cfg.ica.numIC, 'interrupt','off' );  %  ,'stop', 1E-7  );   
    end
    
    % outputs
    cfg.ica.tElapsed = toc(tin_ica); 
    cfg.ica.weights = EEG.icaweights; % store weights reliably in dataset
    cfg.ica.sphere  = EEG.icasphere;
    fprintf('... finished ICA in %.2f mins\n\n\n', cfg.ica.tElapsed/60)
% else
%     fprintf('skipping ICA ...\n')
%     cfg.ica.tElapsed   = nan;
%     cfg.ica.numIC      = nan;
%     cfg.ica.weights    = nan;
%     cfg.ica.sphere     = nan;
% end

%% outputs 
EEG.etc.dusk2dawn.ica  = cfg.ica;
% EEG.etc.dusk2dawn.pica = cfg.pica; 

EEG_out = EEG; % for output of function, not for saving the data to file

% %% save dataset in specified path
% if isfield(cfg, 'savePath')
%     ctemp = []; ctemp.savePath = cfg.savePath; ctemp.headerOnly = true;
%     cleanSleep_saveDataset(EEG, cfg);
% end

end % end function
