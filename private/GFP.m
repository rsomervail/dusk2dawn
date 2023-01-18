%
%
%       compute Global Field Power between 2 eeg samples 
%           - ASSUMES AVG REF
%       'Topographic ERP analyses: A step-by-step tutorial review' - Murray, Brunet & Michel, 2008
%               - Richard Somervail, 2018
%
%
function out = GFP( in )

if max(mean(in)) > 0.1
    error 'RS: data are not average-referenced!';
end

n = length(in);
out = sqrt(  (1/n) * sum(  in.^2   )           );

end