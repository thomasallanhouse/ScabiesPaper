function [chainB, chain_ll, chainTime, ratioA, ratioR, thetaNew] = mcmcScabies081015(k,numRuns)
%
% Function designed to run the Metropolis-Hastings MCMC
%
% Input:
%   h: The tuning parameter. Determines how ambitious I am in the
%   multivariate-normal proposals
%
%   k: The method used to solve the master equation
%   1: Expokit
%   2: Chebyshev expansion
%   3: Runge-Kutta (4,5) i.e. ode45
%   4: DA order 1
%   5: DA order 2
%   6: DA order 3 ... an so on
%
%   numRuns: Length of the MCMC chain
%
% Written by Tim Kinyanjui on 8th Oct 2015
% University of Manchester


% Load the data
data = importdata('nit.txt');

% Get data in the proper vectors
[N, I, T] = deal(data(:,1), data(:,2), data(:,3));

% The starting parameters
beta = 0.006; %0.0047
alpha = 0.4; % 0.663
tau = 0; % No infection from community
gamma = 3.6;
 
% Tuning parameter Order: [beta alpha gamma]
% h = [0.0085 0.2];
h = [0.004 0.2];

% Var-cov matrix calculated from a prelinimary run
covM = [2.01490864662294e-03,0.000923984454369986;0.000923984454369986,0.184034178856011];

% Anonymous function
nLL = @(phat)logScabiesGPU(phat,k,tau,N,I,T);

% Anonymous functions for the priors of beta and alpha respectively
priorB = @(x)log(normpdf(x,0.0047,covM(1,1)));
priorA = @(x)log(normpdf(x,0.663,covM(2,2)));

% Pre-allocate
chain_ll = zeros(numRuns,1);
chainB = zeros(numRuns,2);
chainTime = zeros(numRuns,1);

% Calculate the first likelihood
chain_ll(1) = nLL([beta,alpha,gamma]);
chain_ll(2) = chain_ll(1);
chainB(1,:) = [beta alpha];
chainTime(1) = 1;
thetaOld = [beta alpha];
thetaNew = thetaOld;
ratioR = [0 0];
ratioA = [1 1];

% For all the timepoints, propose and select the newdata
for i = 2:numRuns
    
    % Copy to new
    thetaNew(i,:) = thetaNew(i-1,:);
    
    for j = 1:2
        
        % Store the old likelihood
        tHood = chain_ll(i);
        
        % Propose a new set of variables from normal distribution
        thetaNew(i,j) = chainB(i-1,j) + h(j)*randn;
        
        % Reject if proposed is negative
        if thetaNew(i,j) < 0
            
            % Reject the proposed and keep the old
            thetaNew(i,j) = thetaOld(1,j);
            chain_ll(i) = tHood;
            ratioR(1,j) = ratioR(1,j) + 1;
            
        else
            
            % Prior for beta
            if j == 1
                
                % Prior of current value
                pC = priorB(thetaOld(1,j));
                
                % Prior of proposed value
                pP = priorB(thetaNew(i,j));
                
                % Prior for alpha
            elseif j == 2
                
                % Prior of current value
                pC = priorA(thetaOld(1,j));
                
                % Prior of proposed value
                pP = priorA(thetaNew(i,j));
                
            end
            
            % Likelihood of proposed value
            llhood_new = nLL([thetaNew(i,1),thetaNew(i,2),gamma]);
            
            % Always accept if likelihood is improves
            if rand <= exp(llhood_new + pP - tHood - pC)
                
                % Accept
                thetaNew(i,j) = thetaNew(i,j);
                chain_ll(i) = llhood_new;
                ratioA(1,j) = ratioA(1,j) + 1;
                
            else
                
                % Reject
                thetaNew(i,j) = thetaOld(1,j);
                chain_ll(i) = tHood;
                ratioR(1,j) = ratioR(1,j) + 1;
                
            end
            
        end
        
        % Update the parameters
        thetaOld = thetaNew(i,:);
        chainB(i,:) = thetaNew(i,:); %#ok<*AGROW>
        
    end
    
    % Update time
    chainTime(i) = i;
    chain_ll(i+1) = chain_ll(i);
    
end

% Clip the end. This style isn't good programming to do such things
chain_ll(end) = [];

% End
return