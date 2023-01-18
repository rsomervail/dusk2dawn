%
%
%       apply ASR cleaning using d2d_applyCleaning
% 
%       todo:
%           ! %!!!!! currently will fail if less than 3 pars
%           ! give choice about whether to overwrite raw data in memory
% 
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net
% History:
% 11/07/2022 ver 0.1. Created
% 
%%  

function EEG = pop_d2d_applyCleaning(EEG)

% load varied parameters from dataset
pars = EEG.etc.dusk2dawn.cfg.pars;
npars = length(pars.labels);

% create GUI
geometry = { ...
            1 ... % title
            1 ...
            [1 1] ... % par1 
            1 ...
            [1 1] ... % par2
            1 ...
            [1 1] ... % par3
            1 ...
           };

% DEFINITIONS
space = { {} };

% get parameter value lists if any parameters were varied  
if npars >= 1
    parvals1 = d2d_getParList(pars.values{1});
    row1 = {    { 'Style', 'text', 'string', ['parameter 1: ' pars.labels{1}], 'fontweight', 'bold' } ...
                { 'Style', 'popupmenu', 'string', parvals1, 'value', 1,  'tag' 'sel_par_1'  }; };
else 
    row1 = space;
end
if npars >= 2
    parvals2 = d2d_getParList(pars.values{2});
    row2 = {    { 'Style', 'text', 'string', ['parameter 2: ' pars.labels{2}], 'fontweight', 'bold' } ...
                { 'Style', 'popupmenu', 'string', parvals2, 'value', 1,  'tag' 'sel_par_2'  }; };
else 
    row2 = space;
end
if npars >= 3
    parvals3 = d2d_getParList(pars.values{3});
    row3 = {    { 'Style', 'text', 'string', ['parameter 3: ' pars.labels{3}], 'fontweight_', 'bold' } ...
                { 'Style', 'popupmenu', 'string', parvals3, 'value', 1,  'tag' 'sel_par_3'  }; };
else 
    row3 = space;
end


% run GUI if any varied parameters
if npars > 0

    uilist =[ ...  
            {{ 'Style', 'text', 'string', 'Choose a set of tested parameters to apply to your data:', 'fontweight', 'bold' }} ...
            space ...
            row1 ...
            space ...
            row2 ...
            space ...
            row3 ...
            space ...
            ];

      [ tmp1 tmp2 strhalt cfg ] = inputgui( 'geometry', geometry, 'uilist', uilist, ...
           'helpcom','pophelp(''pop_dusk2dawn'');', 'title', 'Apply ASR cleaning to dataset -- pop_dusk2dawn()');

    % call actual function 
    EEG = d2d_applyCleaning(EEG, cfg);

else
    % call actual function 
    EEG = d2d_applyCleaning(EEG);
end



end