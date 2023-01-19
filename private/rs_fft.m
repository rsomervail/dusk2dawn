% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net
% History:
% 19/01/2023 ver 1.0.0 Created
% 
%%  
%%
%
%      compute FFT for one vector (e.g. a single channel over time)
% 
%       to compute frequencies:
%           fstep = 1/((1/srate)*length(data)); 
%           freqs = 0  : fstep : (length(fft)-1) * fstep; 
%         
%
function [out] = rs_fft(in)

out = abs(fft(in(:,1))); % compute FFT and transform to vector of amplitudes 
out = out/size(out,1); % normalize by N 
% out = out.^2; % compute power
out = out(  1 : fix(size(out,1)/2) ); % only output first half of fft (symmetricaly anyways)


end
