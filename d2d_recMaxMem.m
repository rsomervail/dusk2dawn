%
%
%   subfunction which returns a recommended amount of memory to allocate to the ASR cleaning
% 
%
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net 
%%  
function RecMem = d2d_recMaxMem(useGPU)
    
    try % if can use this function (windows)
        if useGPU
            gpu = gpuDevice;
            RecMem = (gpu.AvailableMemory/1024/1024) * 0.8;
        else
            [~, mem] = memory;  
            RecMem = (mem.PhysicalMemory.Available/1024/1024) * 0.75; % erring on the side of caution to reduce chance of memory crash
        end 
    catch % if can't use memory function (use simple defaults)
        if useGPU
            RecMem = 0.5 * 1024;
        else
            RecMem = 4 * 1024;
        end
    end

end
