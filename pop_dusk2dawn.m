%
%   Usage:
%   >> EEG = pop_dusk2dawn(EEG); % launch a GUI
% 
% ## OVERVIEW ##
% Dusk2Dawn allows you to easily clean whole-night sleep EEG using Artifact 
% Subspace Reconstruction (ASR). As ASR must be calibrated using relatively 
% clean reference data, and there is considerable variability in signal-to-noise 
% ratio across and within sleep-stages, a standard ASR applied to the whole-night 
% would lack sensitivity during lower-amplitude N1 and REM stages, while inappropriately
% removing large-amplitude brain signals such as slow-waves in N2 and N3 stages. 
% 
% Dusk2Dawn therefore implements two methods to solve this problem:      
%       (1) running ASR separately to each sleep stage 
%       (2) running ASR in a short sliding-window (e.g. 4 - 16 minutes).
% 
% The other main benefit of Dusk2Dawn is that it allows you to easily and automatically 
% run ASR several times % with different sets of parameters (e.g. different levels of ASR
% severity) and easily test the effects of each run-through on your data, before choosing 
% which final dataset to use for your analysis. This can be achieved by plotting the data 
% before & after ASR superimposed and also by plotting a range of validation metrics including 
% Slow-Wave amplitude & consistency and reduction of power in particular frequency bands. 
% This functionality is fully implemented in the graphical interface of Dusk2Dawn, so that 
% anyone can perform this validation regardless of programming ability. 
% 
% Both of these features of Dusk2Dawn make it a useful tool to clean whole-night sleep EEG, 
% or indeed any long EEG recordings.
%
% An explanation of how ASR  works can be found on the github page for the clean_rawdata plugin
% for EEGLAB: https://github.com/sccn/clean_rawdata/blob/master/README.md
%
% 
% ## PLEASE CITE THE FOLLOWING PAPER ##
% Somervail R, Cataldi J, Stephan AM, Siclari F, Iannetti GD. 2023. 
% Dusk2Dawn: an EEGLAB plugin for automatic cleaning of whole-night sleep electroencephalogram using Artifact Subspace Reconstruction. 
% Sleep. 1–14.
% 
% 
% ## HOW TO USE ##
% This function calls both D2D core functions and will first clean the data, and then display
% a window which lets you plot the data before and after ASR cleaning, as well as plot several 
% validation metrics, including spectral decomposition of the data, Slow-Wave amplitude & 
% consistency, and ICA quality.
% 
% 
%   ## Handling multiple datasets in a memory efficient way ##
%   All functions in Dusk2Dawn can be run with multiple datasets selected so that you can save 
%   yourself the extra time of manually running everything for each subject.
% 
%   However, since whole-night EEG is very memory intensive, to do this you must change your 
%   EEGLAB preferences so that only one dataset is stored in memory at a time.
% 
%   To handle multiple datasets in a memory efficient way:
% 
%       (1) Go to EEGLAB -> File -> Preferences, and check the box labelled "If set, keep at 
%           most one dataset in memory".
% 
%       (2) Load multiple datasets and then go to EEGlab -> Datasets -> Select Multiple Datasets, 
%           before running Dusk2Dawn cleaning & validation functions.
% 
% ___________________________________________________________________________________________________________
% 
% Step 1 - Clean data with 1 or more sets of ASR parameters & compute validation metrics
% ___________________________________________________________________________________________________________
% 
%   - First you must choose the ASR parameters for the cleaning. The most important parameter here
%     is the cutoff value. This determines the ASR severity, i.e. how much larger artifacts must be 
%     than signal to be removed by ASR, and is expressed in standard deviations (SD). In sleep we 
%     have found values from 30 to 45 SD are appropriate to avoid removing Slow Waves.
%           
%         - You can vary up to 3 parameters in this screen by entering a range of values like this: 
%               [ 30, 35, 40, 45 ]
%               
%         - Advanced settings can also be found by clicking the popup buttons, but for most users 
%           this is unnecessary.
%            
%   - Second, you must choose the ASR implementation. By default the data will be segmented by sleep 
%     stage and ASR performed in a sliding-window. You must specify the events in your dataset which 
%     correspond to sleep stages using the selection box. You can also vary the sliding window length 
%     and overlap, although the default settings are probably fine here.
% 
%   - Third you can choose which additional validation metrics you would like to compute before & after 
%     cleaning. These include:
% 
%         (1) How much variance ASR removed from the data (computed by default)
%                
%         (2) Power by frequency band
%                
%         (3) Slow-wave amplitude and consistency (here you can either select previously-loaded events 
%             in your dataset which correspond to SW events or automatically detect them in each dataset)
%                
%         (4) Quality of Independent Component Analysis (ICA) decomposition using IClabel. This allows you
%             to judge how much ASR has removed artifacts which are difficult to classify using ICA (labelled
%             'unknown' components) and how much of the data is comprised of brain activity.
% ___________________________________________________________________________________________________________
% 
%  Step 2 - Evaluate the results (plot data & validation metrics) and apply final ASR cleaning
% ___________________________________________________________________________________________________________
% 
%   - After the data has been cleaned using ASR with the select set(s) of parameters 
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
%%  

function pop_dusk2dawn % no inputs or outputs because everything evaluated in base

overwrite = '[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, 1:length(EEG)); eeglab redraw;';

% clean data & validate cleaning
cmd_clean = '[EEG, cancelFlag] = pop_dusk2dawn_clean(EEG);';
evalin("base", cmd_clean);
evalin("base", overwrite);

% evaluate results (plot data before/after ASR, plot validation metrics) & apply final cleaning
cmd_evalResults = [ ...
'if ~cancelFlag, ' ...
'EEG = pop_dusk2dawn_evalResults(EEG);' ...
'end'
];
evalin("base", cmd_evalResults);
evalin("base", ['if ~cancelFlag,' overwrite 'end'] );

end
