%
%
%
%   takes varied parameter cell array and returns as a string that can be printed nicely in an inputgui GUI
% 
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net
% History:
% 11/07/2022 ver 0.1. Created
%%
function parlist = d2d_getParList(vals)


parlist = cell(length(vals),1);
for k = 1:length(vals)
    if isnumeric(vals{k})
        if size(vals{k},1) == 1 %  check if there are rows
            parlist{k,1} = num2str(vals{k});
        else % else transpose to prevent strfind errors
            temp = num2str(vals{k}');
            temp = strsplit(temp,' '); % remove excessive whitespaces
            temp = strjoin(temp,', ');
            parlist{k,1} = temp;
        end
    else
        parlist(k,1) = vals(k);
    end
end


% parvals2 = cellfun( @num2str , pars.values{2}, 'UniformOutput',false );
% parvals3 = cellfun( @num2str , pars.values{3}, 'UniformOutput',false );


end