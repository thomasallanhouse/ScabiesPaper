% Functions runs the SEI gillespies model for each of the MCMC sample
% Written by Tim Kinyanjui
% on 8th May 2017

function [Y, fullSet] = SEI_Gillespie_multiple(sim_time,Y0,beta,alpha,gamma)

% Loop for each run
parfor i = 1:length(beta)
    
    % Run the model
    [time, YY] = SEI_Gillespie(sim_time,Y0,beta(i),alpha(i),1/gamma(i));
    
    % Store the important data i.e. sum of E and I at the last time point
    Y(i) = sum(YY(end,2:3));
    
    % Also store the entire matrix
    fullSet{i,1} = [time' YY];

end

% End
return