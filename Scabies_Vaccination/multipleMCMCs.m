% Script runs the MCMC in parallel
clearvars % This is important to ensure no errors
tic;

parfor i = 1:16
    
    % Run the chains
    [chainB(:,:,i), chain_ll(:,i), chainTime(:,i), ratioA(i,:), ratioR(i,:), thetaNew(:,:,i)] = mcmcScabiesAllCombined190917(7,25000);
    
end

timer = toc;

% Save the results
save multipleChains_051017
