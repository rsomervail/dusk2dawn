%
%
%       - takes data split by the d2d_splitByStage function and recombines them into the original continuous data
% 
%       inputs:
%           EEG_all   - a struct array of EEGlab structures containing the separated data for each sleep stage
%           EEG_dummy - a blank version of the original data used as a template for the re-combination of the seperated stages
% 
%               - Richard Somervail, 22/11/2021
%                   www.iannettilab.net
%%
function EEG_out = d2d_recombineByStage( EEG_all , EEG_dummy  )

nstages = length(EEG_all);

% merge validation structures
for st = 1:nstages
    valid(st) = EEG_all(st).etc.dusk2dawn.valid;    
end
for st = 1:nstages
    valid(st).stage = EEG_all(st).etc.dusk2dawn.stageSplit.stageCodes(EEG_all(st).etc.dusk2dawn.stageSplit.thisStage);
end

% merge actual data
EEG_out = EEG_dummy; clear EEG_dummy
for st = 1:nstages
    indy2take = EEG_all(st).etc.dusk2dawn.stageSplit.indy;
    EEG_out.data(:,indy2take) = EEG_all(st).data; 
    EEG_all(st).data = []; % clear this stage 
end

% add validation structure to output
EEG_out.etc.dusk2dawn.valid = valid;



end