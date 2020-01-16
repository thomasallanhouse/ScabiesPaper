function [chainA, chain_ll, chainTime, ratioA, ratioR, thetaNew] = mcmcScabies091015(k,numRuns)
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
beta = 0.0047; %0.01;
alpha = 0.4; % 0.75;
tau = 0; % No infection from community
gamma = 3.6;
h=0.08;

% Anonymous function
nLL = @(phat,red1,red2,red3)logScabiesGPU(phat,k,tau,N,I,T);
prior = @(x)log(unifpdf(x,0,5));

% Calculate the first likelihood
chain_ll(1) = nLL([beta,alpha,gamma]);
chainA(1) = alpha;
chainTime(1) = 1;
thetaOld = alpha;
thetaNew = thetaOld;
ratioR = 0;
ratioA = 1;

% For all the timepoints, propose and select the newdata
for i = 2:numRuns
    
    % Store the old likelihood
    tHood = chain_ll(i-1);
    
    % Propose a new set of variables from normal distribution
    thetaNew(i,1) = chainA(i-1) + h*randn;
    
    % Reject if proposed is negative
    if thetaNew(i,1) < 0
        
        % Reject the proposed and keep the old
        thetaNew(i,:) = thetaOld;
        chain_ll(i) = tHood;
        ratioR = ratioR + 1;
        disp('here')
        
    else
        
        % Prior of current value
        pC = prior(chainA(i-1));
        
        % Prior of proposed value
        pP = prior(thetaNew(i,1));
        
        % Likelihood of proposed value
        llhood_new = nLL([beta,thetaNew(i,1),gamma]);
        
        % Always accept if likelihood is improves
        if rand <= exp(llhood_new + pP - tHood - pC)
            
            % Accept
            thetaNew(i,:) = thetaNew(i,:);
            chain_ll(i) = llhood_new;
            ratioA = ratioA + 1;
            
        else
            
            % Reject
            thetaNew(i,:) = thetaOld;
            chain_ll(i) = tHood;
            ratioR = ratioR + 1;
            
        end
        
    end
    
    % update the parameters
    thetaOld = thetaNew(i,:);
    chainA(i) = thetaNew(i,1); %#ok<*AGROW>
    chainTime(i) = i;
    
end

% End
return