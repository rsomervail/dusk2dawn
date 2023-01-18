%
%
%       !!! WRITE DOCUMENTATION
% 
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net
% History:
% 11/07/2022 ver 0.1. Created
% 
%%  

function EEG = pop_dusk2dawn(EEG)

EEG = pop_dusk2dawn_clean(EEG);

EEG = pop_dusk2dawn_evalResults(EEG);

end