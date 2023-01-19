%
%
%   subfunction which returns a recommended amount of memory to allocate to the ASR cleaning
% 
%
% 
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net
% History:
% 19/01/2023 ver 1.0.0 Created
% 
%%  
function RecMem = d2d_recMaxMem(useGPU)

    if useGPU
        gpu = gpuDevice;
        RecMem = (gpu.AvailableMemory/1024/1024) * 0.8;
    else
        [~, mem] = memory;  
        RecMem = (mem.PhysicalMemory.Available/1024/1024) * 0.75; % erring on the side of caution to reduce chance of memory crash
    end

end
