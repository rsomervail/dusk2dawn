%
%       converts a Nx2 matrix of indices ... where 
%       column 1 = start of each segment of consecutive ones, and column 2
%       = the end of that segment ... into a logical array
%           - i.e. its the inverse of rs_logical2indices
%           * also requires length of original array *
% 
%       out = indices2logical(in, len)
%
%% 

function  out = indices2logical(in, len)

out = false(1,len);
for k = 1:size(in,1)
    out( in(k,1):in(k,2) ) = true;
end % figure; plot(out)


end