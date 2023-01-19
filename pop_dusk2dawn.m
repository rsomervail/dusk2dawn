%
%
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net
% History:
% 19/01/2023 ver 1.0.0 Created
% 
%%  

function EEG = pop_dusk2dawn(EEG)

EEG = pop_dusk2dawn_clean(EEG);

EEG = pop_dusk2dawn_evalResults(EEG);

end
