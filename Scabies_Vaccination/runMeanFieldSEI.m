% Load the data
function [time, Y] = runMeanFieldSEI(t,y0,beta,alpha,gamma)

% Options
options = odeset('RelTol',1e-7);

% Run the model
parfor j = 1:100
    
    [time(:,j), Y(:,:,j)] = ode45(@meanFieldSEI,t,y0,options,beta(j),alpha(j),gamma(j));
    
end

% End
return

% t = [0:1:400];
% y0 = [49, 0, 1];
% load multipleChains_240417
% burnin = 1000;
% thin = 15;
% chainBB = chainB(burnin:thin:end,1,:); beta = chainBB(:);
% chainA = chainB(burnin:thin:end,2,:); alpha = chainA(:);
% chainG = chainB(burnin:thin:end,3,:); gamma = chainG(:);