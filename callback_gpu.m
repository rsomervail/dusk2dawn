% 
% 
% 
%   subfunction of advanced_popup_asr, for updating asr_maxmem paramter in GUI when switching GPU on/off
% 
% 
function callback_gpu
    
    % get new value of asr_useGPU setting
    fig = gcf;
    for k = 1:length(fig.Children)
        if strcmp(fig.Children(k).Tag,'asr_useGPU')
            cur_gpu = fig.Children(k).Value;
            break
        end
    end

    % get recommended or default memory setting
    maxmem = d2d_recMaxMem(cur_gpu);

    % update editbox value
    for k = 1:length(fig.Children)
        if strcmp(fig.Children(k).Tag,'asr_MaxMem')
            fig.Children(k).String = [ num2str(maxmem/1024) '*1024' ];
            break
        end
    end

end