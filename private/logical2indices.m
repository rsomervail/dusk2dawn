%
%       converts a logical array into a Nx2 matrix of indices, where 
%       column 1 = start of each segment of consecutive ones, and column 2
%       = the end of that segment
%
%
%
%
%% 

function  out = logical2indices(in)

%   figure; plot(in)
difs = diff([ 0 in 0 ]);  %   figure; plot(difs)
[ ~, inds1 ] =  find(difs > 0);
[ ~, inds2 ] =  find(difs < 0); 
inds2 = inds2 - 1; % remove one so that this corresponds to the final index of the segment
out = [inds1', inds2'];

end