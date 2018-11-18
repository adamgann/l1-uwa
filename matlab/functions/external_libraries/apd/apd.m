function [a, p] = apd(s);
% [a, p] = apd(s) estimates the amplitude probability distribution function.
%
% input parameters:
% s = amplitude samples
%
% return variables:
% a = ordered amplitudes
% p = probability that the ordered ampitude is exceeded
% R.J. Achatz, et al., US DoC, NTIA/ITS 
if isreal(s) & min(s)>=0
 a = sort(s);
 N = length(a);
 p = 1 - [1:N]/N;
else
 disp('Input values must be amplitudes i.e. real and positive values.');
end 