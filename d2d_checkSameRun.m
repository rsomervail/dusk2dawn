%
% check ASR/D2D settings are the same across all datasets in EEG struct
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net 
function flag = d2d_checkSameRun(EEG)
    for f = 1:length(EEG)
        temp = EEG(f).etc.dusk2dawn.cfg;
        d2d{f} = rmfield( temp, {'origFile','cleanFiles', 'asr_useGPU'} );   
        d2d{f} = rmfield( d2d{f},  );
    end
    flag = isequaln(d2d{:});
    
    if ~flag
        errordlg2('Datasets do not have matching ASR parameters, run again with the same settings ' ...
            ,'Error: datasets cannot be merged');
        error 'Error: datasets are not compatible'
    end
end