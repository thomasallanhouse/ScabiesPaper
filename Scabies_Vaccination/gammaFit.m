function muci = gammaFit(x0)

% Function fits gamma distribution to the data i.e. 95%CI of the mean
% Fit this using lsqcurvefit
% Written by Tim Kinyanjui on 2nd Nov 2015

% Generate data
%data = gamrnd(x0(1),x0(2),10000,1);

% % Fit a normal distribution
% [meanG,varG] = gamstat(x0(1),x0(2));
% 
% muci(1,1) = meanG - (2*sqrt(varG));
% muci(2,1) = meanG + (2*sqrt(varG));

% Generate random data
data = gamrnd(x0(1),x0(2),100000,1);

% Gamfit
phat = gamfit(data);


% Calculate the cdf
lowerP = gamcdf(2,phat(1),phat(2));
upperP = gamcdf(28,phat(1),phat(2));
muci = abs(0.975-(upperP-lowerP));

% End of the function
return

% Run this on the comand window to calculate A and B for the gamma
% distribution
% [xopt, fval] = fminsearch(@gammaFit,[2,5],optimset('Display','iter','TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',5000,'MaxIter',1000))