%
%
%   subfunction which returns a recommended amount of memory to allocate to the ASR cleaning
% 
%               - Richard Somervail, 22/11/2021
%                   www.iannettilab.net
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