%
%       given a raw data filename (origFile)
%           -> get the filename for the cleaned version of data with filename = origFile
%       given cfg with parameter selection fields (sel_par_1, sel_par_2, sel_par_3)
%           -> load the cleaned file with the specified set of varied cleaning parameters 
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022 (function added in Jan 2025)
%           www.iannettilab.net 
%% 
function file2load = d2d_getCleanedFileName(origFile,cfg)

    %% handle inputs
    if ~exist('cfg','var')
        cfg = struct;
    end
    origFile = strrep(origFile,'.set','');

    %% get filename
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
    
end