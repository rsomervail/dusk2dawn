%
% check ASR/D2D settings are the same across all datasets in EEG struct
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net 
function flag = d2d_checkSameRun(EEG)

% concat all d2d settings
    nfiles = length(EEG);
    d2d  = cell(nfiles,1);
    pars = d2d;
    splitByStage = nan(nfiles,1);
    for f = 1:nfiles
        temp = EEG(f).etc.dusk2dawn.cfg;
        d2d{f} = rmfield( temp, {'origFile','cleanFiles', 'asr_useGPU', 'asr_MaxMem', 'savePath'} );  
        splitByStage(f) = d2d{f}.splitByStage;
        pars{f} = temp.pars;
    end
   
    % if split by stage inactive the ignore the stagecodes field
    if sum(splitByStage) ~= nfiles
        for f = 1:nfiles
            d2d{f} = rmfield( d2d{f}, 'stageCodes' );
        end
    end
   
    % check for matching d2d settings
    flag = isequaln(d2d{:});
    if ~flag
        errordlg2('Datasets do not seem to have matching ASR parameters, run again with the same settings ' ...
            ,'Warning: datasets may not be compatible, this may cause errors when plottin');
    end
   
    % check for matching d2d varied parameters (i.e. the most crucial check)
    if ~isequaln(pars{:})
        errordlg2('Datasets do not seem to have matching *varied* ASR parameters, run again with the same settings ' ...
            ,'Error: plotting validation metrics will not be possible if datasets do not share the same varied parameters');
        error 'Error: plotting validation metrics will not be possible if datasets do not share the same varied parameters';
    end

end