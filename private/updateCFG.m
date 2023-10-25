%  
%  
%   Update configuration structure with advanced settings from popup GUI windows.
%   Probably no reason for any user to call this function directly.
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net 
%% 
function cfg = updateCFG(cfg,add)
    flds = fields(add);
    for f = 1:length(flds)
        cfg.(flds{f}) = add.(flds{f});
    end
end