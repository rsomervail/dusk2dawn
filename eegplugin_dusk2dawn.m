% eegplugin_clean_dusk2dawn() - a wrapper to plug-in dusk2dawn into EEGLAB. .
% 
%
%
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net
% 
%%  
function vers = eegplugin_dusk2dawn(fig,try_strings,catch_strings)

vers = '3.4.0';
p = fileparts(which('eegplugin_dusk2dawn'));

% add boundedline subdirectory to matlab path
addpath(genpath([p filesep 'boundedLine']))

% DEFINITIONS
overwrite = '[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); eeglab redraw; ';
% titleCol = [0.5 * ones(1,3)];
titleCol = [0.0, 0, 100]/100; % solid blue
subCol = [49.3, 70.3, 100]/100; % ice blue

% CREATE SUBMENU
menuFolder = findobj(fig, 'tag', 'tools');
submenu = uimenu(menuFolder, 'text', 'Dusk2Dawn - Clean raw data using ASR & validate results of cleaning'...
    ,'separator',0, 'ForegroundColor', titleCol, 'position', 7 ...
    , 'userdata', 'startup:off;epoch:off;study:on' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASTER FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uimenu( submenu, 'text', '- MASTER FUNCTION -', 'Separator', 0, ...
%     'userdata' , 'startup:off;continuous:off;epoch:off;study:off;erpset:off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     pop_dusk2dawn
cmd = [ 'pop_dusk2dawn;' ];
uimenu( submenu, 'text', ...
    'Dusk2Dawn - Run all core functions                       (pop_dusk2dawn)', 'userdata', 'startup:off;epoch:off;study:on', ...
    'callback', cmd, 'Separator', 0, 'ForegroundColor', titleCol);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CORE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uimenu( submenu, 'text', ' ' , 'userdata' , 'startup:off;continuous:off;epoch:off;study:off;erpset:off');
uimenu( submenu, 'text', '- CORE FUNCTIONS -', 'Separator', 1, ...
    'userdata' , 'startup:off;continuous:off;epoch:off;study:off;erpset:off');
ident = '   '; % space before text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) dusk2dawn_clean
cmd = [ '[EEG] = pop_dusk2dawn_clean(EEG);' overwrite ];
uimenu( submenu, 'text', ...
    [ident '(1) Clean whole-night EEG using ASR                  (pop_dusk2dawn_clean)'], 'userdata', 'startup:off;epoch:off;study:on', ...
    'callback', cmd, 'Separator', 1);

% (2) dusk2dawn_evalResults (plotValidation, vis_artifacts and applyCleaning)
cmd = [ '[EEG] = pop_dusk2dawn_evalResults(EEG);' overwrite ];
uimenu( submenu, 'text', ...
    [ident '(2) Evaluate results & apply the cleaning           (pop_dusk2dawn_evalResults)'], 'userdata', 'startup:off;epoch:off;study:on', ...
    'callback', cmd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uimenu( submenu, 'text', ' ' , 'userdata' , 'startup:off;continuous:off;epoch:off;study:off;erpset:off');
uimenu( submenu, 'text', '- SUBFUNCTIONS -', 'Separator', 1,...
    'userdata' , 'startup:off;continuous:off;epoch:off;study:off;erpset:off');
% ident = '    '; % space before text
ident = [ident ident]; % space before text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pop_d2d_applyCleaning 
cmd = ['EEG = pop_d2d_loadData(EEG);' overwrite];
uimenu( submenu, 'text', ...
    [ident '- Apply ASR cleaning / Revert to raw data         (pop_d2d_loadData)'], 'userdata', 'startup:off;epoch:off;study:on', ...
    'callback', cmd, 'Separator', 1);

% pop_d2d_vis_artifacts(EEG) 
uimenu( submenu, 'text', ...
    [ident '- Plot data before/after ASR cleaning               (pop_d2d_vis_artifacts(EEG))'], 'userdata', 'startup:off;epoch:off;study:on', ...
    'callback', 'pop_d2d_vis_artifacts(EEG)');

% pop_d2d_plotValidation   
uimenu( submenu, 'text', ...
    [ident '- Plot effects of ASR on validation metrics        (pop_d2d_plotValidation)'], 'userdata', 'startup:off;epoch:off;study:on', ...
    'callback', 'pop_d2d_plotValidation(EEG)');


% split by stage !!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% final line of free space at the end of the menu
uimenu( submenu, 'text', ' ' , 'userdata' , 'startup:off;continuous:off;epoch:off;study:off;erpset:off');


