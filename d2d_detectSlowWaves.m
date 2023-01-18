%
%       default slow-wave detection function for cleanSleep plugin for EEGlab
%           - uses the "swalldetectnew.m" function and adds the output to the events structure of a EEGlab variable
%           - currently only handles one channel
%           
%               - Richard Somervail, 30/11/2021
% 
function [EEG] = d2d_detectSlowWaves(EEG, chans, ampthreshold)

% run swalldetectnew function
[swa_results] = swalldetectnew( EEG.data(chans,:), EEG.srate, ampthreshold );

% get negative peaks and add to EEG struct events
negpks = swa_results.channels.maxnegpk; 
EEG = eeg_addnewevents(EEG, negpks, repmat({'SW_neg'},length(negpks),1) );

% get positive peaks 
pospks = swa_results.channels.maxpospk; 
EEG = eeg_addnewevents(EEG, pospks, repmat({'SW_pos'},length(pospks),1) );


end