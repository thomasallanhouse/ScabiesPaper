% Functions runs the SEI gillespies model for each of the MCMC sample
% Written by Tim Kinyanjui
% on 8th May 2017

function [Y, fullSet] = SEIR_Gillespie_multiple_Vaccination(sim_time,Y0,beta,alpha,gamma,tVacc,q,vFlag)
% Input:
%   tVacc: Time at whic to treat by a single impulse
%   q:     Probability of a successfull treatment
%   vFlag: 1) if only treating E
%          2) if only treating I
%          3) if treating both E and I
%
% Loop for each run
parfor i = 1:length(beta)
    
    % Rate of loss of infectiousness - Considered impossible in the absence of treatment. Infection not limiting
    tau = 0;
    
    % Run the model
    [time, YY] = SEIR_Gillespie_Vaccination(sim_time,Y0,beta(i),alpha(i),1/gamma(i),tau,tVacc,q,vFlag);
    
    % Store the important data i.e. sum of E and I at the last time point
    Y(i) = sum(YY(end,2:3));
    
    % Also store the entire matrix
    fullSet{i,1} = [time YY];

end

% End
return