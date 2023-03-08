%  
%  
%   Gets the setings used to clean another dataset and applies them to the selected datasets
%     (called from the pop_dusk2dawn_clean GUI
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net 
%%  
function d2d_getSettingsFromSet

    
    % choose file & get cfg from the chosen file
    [filename, filepath] = uigetfile('MultiSelect','off','*.set');
    info = pop_loadset('filename',filename,'filepath',filepath,'loadmode','info');
    cfg = info.etc.dusk2dawn.cfg;
    
    % remove cfg options that relate to the loaded file
    cfg = rmfield( cfg, {'origFile','cleanFiles'} );
    cfg.stageCodes(strcmp(cfg.stageCodes,'[unscored]')) = [];

    % loop through CFG options and update figure accordingly ? could probs condense literally any option here to an if statement on "style"
    fig = gcf;
    for k = 1:length(fig.Children) % ? I could probably manually update all options inside this figure here, not sure about popups though
%         fig.Children(k).Tag % for plotting fields
        switch fig.Children(k).Tag

            case 'splitBySlidingWindow'
                fig.Children(k).Value = cfg.splitBySlidingWindow;
            case 'splitByStage'
                fig.Children(k).Value = cfg.splitByStage;
            case 'asr_cutoff'
                fig.Children(k).String = [ '[ ' num2str(cfg.asr_cutoff) ' ]' ];
            case 'stageWin'
                fig.Children(k).String = [ '[ ' num2str(cfg.stageWin) ' ]' ];
            case 'stageCodes'
                % ? would have to compare strings of unique({info.event.type}) codes and find the correct ones from cfg.stageCodes
            case 'fft_run'
                fig.Children(k).Value = cfg.fft.run;
            case 'sw_run'
                fig.Children(k).Value = cfg.sw.run;
            case 'sw_autoFind'
                fig.Children(k).Value = cfg.sw.autoFind;
            case 'sw_codes' % ? would have to compare strings of unique({info.event.type}) codes and find the correct ones from cfg.sw.codes
%                 if isempty(cfg.sw.codes)
%                     fig.Children(k).Value = []; 
%                 else
% 
%                 end
            case 'ica_run'
                fig.Children(k).Value = cfg.ica.run;
            case 'chunk_len'
                fig.Children(k).String = [ '[ ' num2str(cfg.chunk_len) ' ]' ];
            case 'savePath'
                fig.Children(k).String = cfg.savePath;
%             otherwise 
%                 switch fig.Children(k).Style
%                     case 'edit'
% 
%                     case 'check'
%                         
%                 end
        end
    end

    % handle also advanced settings - ASR
    evalin("caller", sprintf(' cfg_asr.ref_tolerances = ''%s''    ;', [ '[ ' num2str(cfg.ref_tolerances) ' ]' ])  );
    evalin("caller", sprintf(' cfg_asr.ref_maxbadchannels = ''%s''    ;', [ '[ ' num2str(cfg.ref_maxbadchannels*100) ' ]' ])  );
    evalin("caller", sprintf(' cfg_asr.asr_windowlength = ''%s''    ;', [ '[ ' num2str(cfg.asr_windowlength) ' ]' ])  );
    evalin("caller", sprintf(' cfg_asr.asr_useGPU = %s;', num2str(cfg.asr_useGPU) )  );
    evalin("caller", sprintf(' cfg_asr.asr_MaxMem = ''%s''    ;', [ '[ ' num2str(cfg.asr_MaxMem) ' ]' ])  );
    evalin("caller", sprintf(' cfg_asr.chunk_overlap = ''%s''    ;', [ '[ ' num2str(cfg.chunk_overlap) ' ]' ])  );
    
    % handle also advanced settings - valid
    % bin freqs
    for k = 1:length(cfg.fft.binFreqs)
        binFreqs{k} = [num2str(cfg.fft.binFreqs(k,1)) ',' num2str(cfg.fft.binFreqs(k,2))];
    end
    binFreqs = [ '[ ' strjoin(binFreqs,'; ') ' ]' ];
    evalin("caller", sprintf(' cfg_valid.fft_binFreqs = ''%s''    ;', binFreqs) );
    % bin freq labels
    binFreqsLabels = [ '{ ''''' strjoin(cfg.fft.binFreqsLabels,'''''; ''''') ''''' }' ];
    evalin("caller", sprintf(' cfg_valid.fft_binFreqsLabels = ''%s''    ;', binFreqsLabels) );
    % sw - peak window
    evalin("caller", sprintf(' cfg_valid.sw_peakwin = ''%s''    ;', [ '[ ' num2str(cfg.sw.peakwin*1000) ' ]' ])  );
    % sw - channel for sw detection
    evalin("caller", sprintf(' cfg_valid.sw_chan = %d  ;', cfg.sw.chan ) );
    % sw - amplitude threshold
    evalin("caller", sprintf(' cfg_valid.sw_ampThresh = ''%d'' ;', cfg.sw.ampThresh) );
    % ica - number of ICs
    evalin("caller", sprintf(' cfg_valid.ica_numIC = ''%d'' ;', cfg.ica.numIC) );


end % function
