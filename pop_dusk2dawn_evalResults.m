%
%   Usage:
%     >> EEG = pop_dusk2dawn_evalResults(EEG);
%
%   - After the data have been cleaned using ASR with the select set(s) of parameters 
%     (using 'pop_dusk2dawn_clean'), you can visualise the data before & after ASR, 
%     and plot the computed validation metrics.
% 
%   - Once you are happy with your ASR parameters this function will then either:
%      -> load the clean data with the selected ASR parameters (click "OK") 
%      -> or revert back to the raw data (click "CANCEL").
%
%
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net
% 
%%  

function EEG = pop_dusk2dawn_evalResults(EEG)

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

if nfiles > 1
    d2d_checkSameRun(EEG);
end

cfg = EEG(1).etc.dusk2dawn.cfg;
pars = cfg.pars; npars = length(pars.labels);

%%% Define GUI elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% title / help text
if npars == 1, plural = ''; else, plural = 's'; end
if npars == 0, npars_str = 'no'; else, npars_str = num2str(npars); end
if nfiles == 1
    titleText{1} = [ 'This dataset has been cleaned using Dusk2Dawn ' sprintf( 'with %s varied parameter%s.', npars_str, plural)];
else
    titleText{1} = ['These datasets have been cleaned using Dusk2Dawn ' sprintf( 'with %s varied parameter%s.', npars_str, plural)];
end
    titleText{2} = 'You can now explore the effects of ASR on your data using the functions below.';
if npars > 0
    titleText{3} = 'Finally, click "OK" to apply ASR with the selected parameters';
else
    titleText{3} = 'Finally, click "OK" to apply ASR ';
end
titleText{4} = 'or click "CANCEL" to exit without applying the ASR (you can always return to this screen later)';
geo_titleText = { ...
    1 ...
    1 ...
    1 ...
    1 ...
    1 ...
};
ui_titleText = { ...
    { 'Style', 'text', 'string', titleText{1}, 'fontweight', 'bold' } ...
    {} ...
    { 'Style', 'text', 'string', titleText{2}, 'fontweight', 'bold' } ...
    { 'Style', 'text', 'string', titleText{3}, 'fontweight', 'bold' } ...
    { 'Style', 'text', 'string', titleText{4}, 'fontweight', 'bold' } ...
}; 

%% function call buttons
geo_funCalls = { ...
    [1 3] ...
    1 ...
    [1 3] ...
};
ui_funCalls = { ...
    { 'Style','pushbutton', 'string','Plot Data', 'callback','pop_d2d_vis_artifacts(EEG)' } ... 
        { 'Style', 'text', 'string', ' --> Visualise data to compare before/after ASR using vis_artifacts() function' } ...
    {} ...
    { 'Style','pushbutton', 'string','Plot Validation Metrics', 'callback','pop_d2d_plotValidation(EEG)' } ...
        { 'Style', 'text', 'string', ' --> Plot various metrics which quantify the effects of ASR on your data' } ...
};

%% runInfo (i.e. info from CFG about each parameter/setting, whether varied or constant)
%!!! title in italics, everything slightly indented, 2 columns: parameter and value (can be 'varied')



%% Define GUI elements - choose pars to apply ASR with if any were varied
if npars > 0
    ui_choosePars(1) = {{ 'Style', 'text', 'string', 'Choose which set of ASR parameters to apply to your data:', 'fontweight', 'bold' }};
    geo_choosePars = {1};
    for p = 1:3
        if npars >= p
            parvals = d2d_getParList(pars.values{p});
            ui_choosePars = [  ui_choosePars ...
                            {{ 'Style', 'text', 'string', ['parameter ' num2str(p) ': ' pars.labels{p}], 'fontweight', 'bold', 'HorizontalAlignment','right' }} ...
                                {{ 'Style', 'popupmenu', 'string', parvals, 'value', 1,  'tag' ['sel_par_' num2str(p)]  }} ...
                                { {} } ...
                            { {} } ...
            ];
            geo_choosePars = [ geo_choosePars, {[1 1 1]}, {1} ];
        end
    end
else
    geo_choosePars = [];
    ui_choosePars  = [];
end

%% build GUI
space = { {} };

% assemble final UI list & geometry
geometry = [       ...
    1 ...
    geo_titleText  ...
    1 ...
    geo_funCalls   ...
    1 ...
    geo_choosePars ...
    1 ...
];
uilist = [         ...
    space          ...
    ui_titleText   ...
    space          ...
    ui_funCalls    ...
    space          ...
    ui_choosePars  ...
    space          ...
];

%% create GUI

[ ~,~, strhalt, cfg ] = inputgui( 'geometry', geometry, 'uilist',uilist, ...
    'title', 'Evaluate results & apply cleaning -- pop_dusk2dawn_evalResults()', ...
   'helpcom','pophelp(''pop_dusk2dawn_evalResults'');' );

 
%% either apply ASR or exit without applying depending on whether click OK or CANCEL
if strcmp(strhalt,'retuninginputui')
    for f = 1:nfiles
        EEG(f) = d2d_loadData(EEG(f),cfg);
    end
end

end
