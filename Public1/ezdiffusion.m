%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab implementation of the EZ-diffusion model
% Reference: Wagenmakers et al. (2007). An EZ-diffusion model for response
% time and accuracy. Psychonomic Bulletin & Review, 14(1), 3-22.
%
% This code estimates the parameters of the EZ-diffusion model from the
% following inputs:
%   - Pc: Proportion correct
%   - VRT: Response time variance
%   - MRT: Mean response time
%   - n: Number of trials (for edge correction)
%
% The outputs are the model parameters:
%   - v: Drift rate
%   - a: Boundary separation
%   - Ter: Nondecision time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [v, a, Ter] = ezdiffusion(Pc,VRT,MRT,n)
s = 0.1;
s2 = s^2;  % The default value for the scaling parameter s equals 0.1

% If Pc equals 0, 0.5, or 1, the method will not work, and an edge
% correction is required
if Pc == 0
    error('Oops, Pc==0')
elseif Pc == 0.5
    disp('Oops, Pc==0.5')
    Pc = 0.5-(0.5/(2*n));  % Edge correction
elseif Pc == 1
    disp('Oops, Pc==1')
    Pc = 1-(1/(2*n));  % Edge correction
end

L = log(Pc/(1-Pc));  % Calculate the logit
x = L*(L*Pc^2-L*Pc+Pc-.5)/VRT;
v = sign(Pc-.5)*s*x^(1/4);  % drift rate
a = s2*L/v;  % boundary separation
y = -v*a/s2;
MDT = (a/(2*v))*(1-exp(y))/(1+exp(y));
Ter = MRT-MDT;  % nondecision time
end
