function [chainB, chain_ll, chainTime, ratioA, ratioR, thetaNew] = mcmcScabiesMultiVar(k,numRuns)
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
%load('optimalFmincon','x','hess','grad')

% Get data in the proper vectors
[N, I, T] = deal(data(:,1), data(:,2), data(:,3));

% The starting parameters
beta = 0.006; %0.0047
alpha = 0.4; % 0.663
tau = 0; % No infection from community
gamma = 3.6;
% h = [0.0015 0.08]; % Order: [beta alpha]
h = [0.009 0.2];
%coVar = [1.14616267357871e-05,0.000333002840114450;0.000333002840114450,0.0134705568594100];
coVar = [1.78221097773056e-05,0.000626368667849704;0.000626368667849704,0.0281553429817819];

% Anonymous function
nLL = @(phat)logScabiesGPU(phat,k,tau,N,I,T);

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
ratioR = 0;
ratioA = 1;

% For all the timepoints, propose and select the newdata
for i = 2:numRuns
    
    % Propose a new set of variables from bi-variate normal distribution
    thetaNew(i,:) = mvnrnd(chainB(i-1,:),coVar);
    
    % Reject if proposed is negative
    while sum(thetaNew(i,:) < 0) > 0
        
        % Propose a new set of variables from bi-variate normal distribution
        thetaNew(i,:) = mvnrnd(chainB(i-1,:),coVar);
        
    end
    
    % Likelihood of proposed value
    llhood_new = nLL([thetaNew(i,1),thetaNew(i,2),gamma]);
    
    % Always accept if likelihood is improves
    if rand <= exp(llhood_new - chain_ll(i-1))
        
        % Accept
        chainB(i,:) = thetaNew(i,:);
        chain_ll(i) = llhood_new;
        ratioA(1) = ratioA(1) + 1;
        
    else
        
        % Reject
        chainB(i,:) = thetaOld(1,:);
        chain_ll(i) = chain_ll(i-1);
        ratioR(1) = ratioR(1) + 1;
        
    end
    
    % Update the parameters and time
    thetaOld = thetaNew(i,:);
    chainTime(i) = i;
    
end

% Clip the end: Not very good programming to do such things
chain_ll(end) = [];

% End
return