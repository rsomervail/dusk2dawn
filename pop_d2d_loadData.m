%
%
%       load a version of the current dataset that has been processed with Dusk2Dawn
%           - e.g. apply ASR cleaning with a particular set of parameters
%              or revert back to the raw data
%
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net 
%%  

function EEG = pop_d2d_loadData(EEG)

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

% load varied parameters from dataset
pars = EEG(1).etc.dusk2dawn.cfg.pars;
npars = length(pars.labels);

%% if no varied parameters then toggle between clean data and raw data
if npars == 0
    cfg = [];
    cfg.loadRaw = ~EEG(1).etc.dusk2dawn.valid(1).raw;
    for f = 1:nfiles

        % toggle between raw and clean data % if multiple sets then keep them all in line
        EEG(f) = d2d_loadData(EEG(f), cfg);
    end
end

%% define GUI elements - choose ASR parameters
if npars > 0
    ui_choosePars(1) = {{ 'Style', 'text', 'string', 'Choose which ASR parameters to apply to your data:', 'fontweight', 'bold' }};
    geo_choosePars = {1};
    for p = 1:3
        if npars >= p
            parvals = d2d_getParList(pars.values{p});
            ui_choosePars = [  ui_choosePars ...
                            {{ 'Style', 'text', 'string', ['parameter ' num2str(p) ': ' pars.labels{p}], 'fontweight', 'bold', 'HorizontalAlignment','right' }} ...
                                {{ 'Style', 'popupmenu', 'string', parvals, 'value', 1,  'tag' ['sel_par_' num2str(p)]  }} ...
                            { {} } ...
            ];
            geo_choosePars = [ geo_choosePars, {[1 1]}, {1} ];
        end
    end
else
    geo_choosePars = [];
    ui_choosePars  = [];
end

%% define GUI elements - or choose raw data checkbox

if npars > 0 
    ui_orRaw = [ ...
        { { 'Style', 'checkbox', 'string' 'OR check this box to revert back to raw data' 'value',0, 'tag' 'loadRaw' } } ...    
        { {} } ...
    ];
    geo_orRaw = {1,1};
end

%% run GUI if any varied parameters
space = { {} };

if npars > 0

    % build GUI
    geometry = [       ...
        geo_choosePars ...
        1 ...
        geo_orRaw      ...
        1 ...
    ];
    uilist = [         ...
        ui_choosePars  ...
        space          ...
        ui_orRaw       ...
        space          ...
    ];

    % run GUI
      [ tmp1 tmp2 strhalt cfg ] = inputgui( 'geometry', geometry, 'uilist', uilist, ...
           'helpcom','pophelp(''pop_dusk2dawn'');', 'title', 'Apply ASR cleaning to dataset -- pop_dusk2dawn()');

    % return if cancelled
    if isempty(cfg)
        return
    end

    % call actual function 
    for f = 1:nfiles
        if isempty(EEG(f).etc.dusk2dawn.cfg.savePath)
            cfg.keepList = true; % if running in memory, then will need to keep the list of EEG data structures
        end
        EEG(f) = d2d_loadData(EEG(f), cfg);
    end
end



end
