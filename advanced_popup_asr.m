%
%   popup menu for advanced asr options in dusk2dawn_clean 
%       - takes input cfg which might be empty or might have existing advanced settings from already opening the GUI
%       - if cancel, then doesn't save settings
%       - if ok, returns settings back to pop_dusk2dawn_clean
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net 
function   out = advanced_popup_asr(in, plot, misc)

    %% if plot is selected then plot the actual GUI, else just return default parameters
    if plot

        % definitions
        indent1 = 0.05;
        veditbox = 1.2;
        v1 = 0.5;

        % build gui  ? will probably have to split this up the way I do for the top function for clarity
        geolocal = { ...
            1 ...  % main title
            1 ... 
            1 ... % ref tolerances title
            [ indent1 0.5 1 ] ... % ref tolerances editbox + help 1
            [ indent1 0.5 1 ] ... % help 2
            1 ... % max bad chans title
            [ indent1 0.5 1 ] ... % max bad chans editbox + help 1
            [ indent1 0.5 1 ] ... % help 2
            1 ... % asr_windowlength title
            [ indent1 0.5 1 ] ... % asr_windowlength editbox + help 1
            [ indent1 0.5 1 ] ... % help 2
            1 ... 
            [ indent1 0.5 1 ] ... % asr_useGPU checkbox + help 1
            1 ... 
            1 ... % advanced d2d options title
            1 ... % sliding window overlap title
            [ indent1 0.5 1 ] ... % sliding window overlap edit box + help 1
            [ indent1 0.5 1 ] ... % help 2
        };
        vertlocal = [ ...
            1 ...  % main title
            v1 ... 
            1 ... % ref tolerances title
            veditbox ... % ref tolerances editbox + help 1
            1 ... % ref tolerances help 2
            1 ... % max bad chans title
            veditbox ... % max bad chans editbox + help 1
            1 ... % help 2
            1 ... % asr_windowlength title
            veditbox ... % asr_windowlength editbox + help 1
            1 ... % help 2
            1 ... % 
            1 ... % useGPU
            1 ... %
            1 ... % advanced d2d options title
            1 ... % sliding window overlap title 
            veditbox ... % sliding window overlap edit box + help 1
            1 ... % help 2
        ];
        uilocal = { ...
            { 'Style', 'text', 'string', 'Advanced ASR Options','fontweight','bold'} ...
            {} ...
            { 'Style', 'text', 'string', 'ref_tolerances - tolerances for acceptable calibration data: [min,max] (SD of RMS signal)'} ...
            {} ... 
                { 'Style', 'edit', 'string', in.ref_tolerances 'tag' 'ref_tolerances'  } ...
                { 'Style', 'text', 'string', '? ASR will automatically find data for calibration whose root-mean-square','fontangle','italic'} ...
            {} ...
                {} ...
                { 'Style', 'text', 'string', '   power does not exceed these tolerance values (e.g. by 5 SD)','fontangle','italic'} ...
            { 'Style', 'text', 'string', 'ref_maxbadchannels - % max bad channels acceptable for calibration data:'} ...
            {} ...
                { 'Style', 'edit', 'string', in.ref_maxbadchannels 'tag' 'ref_maxbadchannels'  } ... 
                { 'Style', 'text', 'string', '? what percentage of channels can exceed the tolerance values above before','fontangle','italic'} ...
            {} ...
                {} ...
                { 'Style', 'text', 'string', '   that segment of data is excluded from the calibration','fontangle','italic'} ...
            { 'Style', 'text', 'string', 'asr_windowlength - length of the sliding-window used by the ASR algorithm when cleaning (s)'} ...
            {} ...
                { 'Style', 'edit', 'string', in.asr_windowlength 'tag' 'asr_windowlength'  } ... 
                { 'Style', 'text', 'string', '? Should not be much longer than the timescale over which artifacts persist, min = 1.5 x numchans / srate','fontangle','italic'} ...
            {} ...
                {} ...
                { 'Style', 'text', 'string', '   (note: this is not the optional D2D sliding window in which both cleaning and calibration are performed)','fontangle','italic'} ...
            {} ...
            {} ...
                { 'Style', 'checkbox', 'string', 'use GPU for ASR processing?' 'value' in.asr_useGPU 'tag' 'asr_useGPU'  } ... 
                { 'Style', 'text', 'string', '? will usually speed up ASR processing, but requires an appropriate graphics card','fontangle','italic'} ...
            {} ...
            { 'Style', 'text', 'string', 'Advanced Dusk2Dawn Options','fontweight','bold'} ...
            { 'Style', 'text', 'string', 'choose overlap of sliding windows (mins): '} ...
            {} ...
                { 'Style', 'edit', 'string', in.chunk_overlap 'tag' 'chunk_overlap'  } ... 
                { 'Style', 'text', 'string', '? something like 1 minute is probably fine, but you can vary this if unsure','fontangle','italic'} ...
            {} ...
                {} ... 
                { 'Style', 'text', 'string', '   although note that larger values result in longer computation time','fontangle','italic'} ...
        }; 

        %% run gui 
        [ ~, ~, strhalt_local, cfg_local ] = inputgui( 'geometry',geolocal, 'uilist',uilocal, 'geomvert',vertlocal, ...
           'title','Advanced Options - ASR Cleaning',  'helpcom', 'pophelp(''pop_dusk2dawn'');'); 
    
        % if clicked OK then return new settings
        if strcmp(strhalt_local,'retuninginputui') 
            out = cfg_local;
        else
            out = in; % else don't change existing settings (return the same)
        end

    %% don't plot just return defaults 
    else
        % ASR defaults
        out.ref_tolerances = '[  -3.5,  5  ]';
        out.ref_maxbadchannels = '[   7.5   ]';
        out.asr_windowlength = ['[ '  num2str( max( 2,1.5*misc.nbchan/misc.srate) , 8 )  ' ]'];
        out.asr_useGPU = 0;

        % d2d defaults
        out.chunk_overlap = '[ 1 ]';
    end

end 
