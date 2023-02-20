%
%   GUI for the wrapper function for the vis_artifacts function in the clean_rawdata plugin for EEGLAB
%
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net 
%%  
function pop_d2d_vis_artifacts(EEG)

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

cfg = EEG(1).etc.dusk2dawn.cfg;
pars = cfg.pars; npars = length(pars.labels);

%% predefinitions
space = { {} };

%% define GUI elements - choose which dataset to plot
if nfiles > 1
    geo_whichData = [ ...
        1 ...
        1 ...
    ];
    ui_whichData = { ...
        { 'Style', 'text', 'string', 'Choose which dataset to plot:', 'fontweight', 'bold' } ...
        { 'Style', 'popupmenu', 'string', {EEG.setname}, 'value', 1,  'tag' 'whichData'  } ...
    };
else
    geo_whichData = [];
    ui_whichData  = [];
end


%% define GUI elements - choose parameters of datasets to plot

if npars > 0

    for p = 1:3
        if npars >= p
            parvals = d2d_getParList(pars.values{p});
            parlabel = { 'Style', 'text', 'string', ['parameter ' num2str(p) ': ' pars.labels{p}], 'fontweight', 'normal'};
            choosePars{p} = {  % row1
                        parlabel ...
                        parlabel ...
                        ... % row2
                        { 'Style', 'popupmenu', 'string', [{' '}; parvals], 'value', 1,  'tag' ['selR_par_' num2str(p)]  } ...
                            {} ...
                            {} ...
                            { 'Style', 'popupmenu', 'string', parvals, 'value', 1,  'tag' ['selL_par' num2str(p)]   } ...
                            {} ...
                            {} ...
            };
        else 
            choosePars{p} = [ repmat(space,1,8)  ];
        end  
    end
    geo_choosePars = { ...
        1 ... % title
        1 ... % title
        1 ...
        [1 1] ... % left/right red/blue
        [1 1] ... % raw data checkbox
        1 ...
        [1 1 ] ... % param 1 (top)
        ones(1,6) ... % param 1 (bottom)
        1 ... % 
        [1 1 ] ... % param 2 (top)
        ones(1,6)  ... % param 2 (bottom)
        1 ... % 
        [1 1 ] ... % param 3 (top)
        ones(1,6)  ... % param 3 (bottom)
    };
    ui_choosePars = [ ...  
        {{ 'Style', 'text', 'string', 'Choose which versions of the dataset to plot together:', 'fontweight', 'bold' }} ...
        {{ 'Style', 'text', 'string', '-> when you are done click "OK"', 'fontweight', 'bold' }} ...
        space ...
        {{ 'Style', 'text', 'string', 'plot in red', 'fontangle','italic' ,'fontweight', 'bold' }} ...
            {{ 'Style', 'text', 'string', 'plot in blue', 'fontangle','italic' ,'fontweight', 'bold' }} ...
        {{ 'Style', 'checkbox', 'string' 'select raw data' 'value' 1 'tag' 'loadRaw' }} ...
            space ...
        space ...
        choosePars{1}  ...
        space ...
        choosePars{2}  ...
        space ...
        choosePars{3}  ...
    ];
else % no parameters to vary
    geo_choosePars = [];
    ui_choosePars  = [];
end


%% build GUI

% assemble final UI list & geometry
geometry = [       ...
    1 ...
    geo_whichData  ...
    1 ...
    geo_choosePars ...
    1 ...
];
uilist = [         ...
    space          ...
    ui_whichData   ...
    space          ...
    ui_choosePars  ...
    space          ...
];

%% create GUI
if (npars > 0) || (nfiles > 1) % these are the two conditions in which we need to plot the GUI here

    [ tmp1 tmp2 strhalt cfg ] = inputgui( geometry, uilist, ...
       ['pophelp(''' mfilename ''');'], ['Superimpose the same dataset before/after ASR cleaning -- ' mfilename '()']);
     
    % return if cancelled
    if ~strcmp(strhalt,'retuninginputui') 
        return
    end

    %% GUI outputs
    % subtract the initial blank space from each selection
    flds = fields(cfg);
    for k = 1:length(flds)
        if startsWith( flds{k}, 'sel_r')
            cfg.(flds{k}) =  cfg.(flds{k}) - 1;
        end
    end

end

%% call d2d_vis_artifacts
if nfiles > 1
    EEG = EEG(cfg.whichData);
end
EEG.data = []; % no need to send the data to lower level functions (it's loaded there anyway)
d2d_vis_artifacts(EEG, cfg);

end
