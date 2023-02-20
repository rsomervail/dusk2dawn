%
%   GUI for d2d_plotValidation & d2d_group_plotValidation which plots the previously computed validation metrics across ASR parameter sets (and possibly datasets)
%
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net 
%%  

function pop_d2d_plotValidation(EEG)

%% get info about EEG data
nfiles = length(EEG);
if any( arrayfun( @(x) ~isfield(x.etc,'dusk2dawn'), EEG ) )
    if nfiles == 1
        errordlg2('Please use the master function (pop_dusk2dawn) or core function (1) first to clean the dataset before running this function' ...
            ,'Error: Dataset has not been cleaned with Dusk2dawn');
    else
        errordlg2('Please use the master function (pop_dusk2dawn) or core function (1) first to clean all the datasets before running this function' ...
            ,'Error: One or more datasets have not been cleaned with Dusk2dawn');
    end
    return
end

% check datasets were cleaned using the same parameters and can therefore be plotted together
if nfiles > 1
    d2d_checkSameRun(EEG);
end

cfg = EEG(1).etc.dusk2dawn.cfg; 
pars = cfg.pars; npars = length(pars.labels);

%!!! DO SAME THING AS VIS ARTIFACTS FUNCTION HERE, JUST REMOVE THE BITS FOR DOING VALIDATION METRIC RANGES IF NO VARIED PARAMS

%% Define GUI elements - selection of stages to plot
if cfg.splitByStage

    stageCodes = [EEG(1).etc.dusk2dawn.valid_merged(:).stage];
    nstages = length(stageCodes);

    geo_stageSel = { ...
        1 ... % title
        1 ...
    };
    vert_stageSel = [ ...
        1 ...
        nstages ...
    ];
    ui_stageSel = {  ...
        { 'Style','text', 'string','choose which sleep stages to plot metrics for', 'fontweight', 'bold' } ...
        { 'Style','listbox', 'string',stageCodes, 'min',1,'max',nstages,'value',1:nstages 'tag','stages' } ... 
    };

else
    geo_stageSel  = {};
    vert_stageSel = [];
    ui_stageSel   = {};
end

%% Define GUI elements - selection of parameters to plot
labels = [{''}, pars.labels];

% determine default view
if npars >= 1, xaxis_default = 2; else, xaxis_default = 1; end
if npars >= 2, raxis_default = 3; else, raxis_default = 1; end
if npars >= 3, caxis_default = 4; else, caxis_default = 1; end

geo_parSel = { ...
    1 ... % title
    1 ...
    [1 1 1]    ... % subtitles
    [1 1 1]    ... % popup menus
};
vert_parSel = [
    1 ...
    1 ...
    1 ...
    1 ...
];
ui_parSel = {  ...
    { 'Style', 'text', 'string', 'choose which parameters to plot on each axis', 'fontweight', 'bold' } ...
    {} ...
    { 'Style', 'text', 'string', 'plot on x-axis', 'fontweight', 'bold' } ...
        { 'Style', 'text', 'string', 'plot as rows', 'fontweight', 'bold' } ...
        { 'Style', 'text', 'string', 'plot as columns', 'fontweight', 'bold' } ...
    { 'Style', 'popupmenu', 'string', labels, 'value', xaxis_default,  'tag' 'plotX'  } ...
        { 'Style', 'popupmenu', 'string', labels, 'value', raxis_default,  'tag' 'plotRows'  } ...
        { 'Style', 'popupmenu', 'string', labels, 'value', caxis_default,  'tag' 'plotCols'  } ...
};

%% Define GUI elements - checkboxes

geo_checkBoxes = [ ...
    1 ...
];
ui_checkBoxes = { ...
    { 'Style', 'checkbox', 'string' ' save figures after plotting?' 'fontweight', 'bold' 'value' 0 'tag' 'saveFig' } ...
};

%% Build GUI
space = { {} };

% assemble final UI list & geometry
geometry = [    ...
    1           ...
    geo_parSel  ...
    1           ...
    geo_stageSel ...
    1           ...
    geo_checkBoxes
];
vert = [ ...
    1 ...
    vert_parSel ...
    1 ...
    vert_stageSel ...
    1 ...
    ones(1,length(ui_checkBoxes));
];
uilist = [         ...
    space          ...
    ui_parSel      ...
    space          ...
    ui_stageSel    ...
    space          ...
    ui_checkBoxes  ...
];

%% create GUI
[ tmp1 tmp2 strhalt cfg ] = inputgui( 'geometry',geometry,'uilist',uilist, 'geomvert',vert, ...
   'helpcom','pophelp(''pop_d2d_group_plotValidation'');','title','Plot effects of ASR on various validation metrics -- pop_d2d_group_plotValidation()');
 
% return if cancelled
if ~strcmp(strhalt,'retuninginputui') 
    return
end

%% extract inputs to d2d_plotValidation from cfg

% handle varied parameters
if isfield(cfg,'plotX')
    cfg.plotX = cfg.plotX - 1;
else
    cfg.plotX = 0;
end
if isfield(cfg,'plotRows')
    cfg.plotRows = cfg.plotRows - 1;
else
    cfg.plotRows = 0;
end
if isfield(cfg,'plotCols')
    cfg.plotCols = cfg.plotCols - 1;
else
    cfg.plotCols = 0;
end
% arrange dimensions to plot
cfg.dims2plot = [ cfg.plotX cfg.plotRows cfg.plotCols  ];
cfg.dims2plot = cfg.dims2plot( [cfg.dims2plot]~=0 ); % remove 0 dims (which are not plotted)
% confirm no parameters have been repeated
if length(cfg.dims2plot) ~= length(unique(cfg.dims2plot)) 
    errordlg2('You have selected the same varied parameter to be plotted in multiple axes - retry and only select 1 axis for each varied parameter',...
        'Error: pop_d2d_plotValidation: Cannot plot same varied parameter on multiple dimensions')
end


%% call subfunction
d2d_plotValidation(EEG, cfg);

end
