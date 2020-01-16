% Plot the Bifurcation diagram
% Written by Tim Kinyanjui
% on 3rd May 2017
% University of Manchester

% Clear workspace
clearvars;

% Determine the parameters to use from the MCMC results
%load multipleChains_081217
load multipleChains_240417
burnin = 1000;
thin = 15;
chainBB = chainB(burnin:thin:end,1,:); beta = chainBB(:);
chainA = chainB(burnin:thin:end,2,:); alpha = chainA(:);
chainG = chainB(burnin:thin:end,3,:); gamma = chainG(:);

% Size of the care homes
N = 1:1:100;

% Initialize the matrices
EI = zeros(length(N),length(N));
initialI = EI;

for j = 1:length(N)
    
    for k = 1:N(j)
        
        % Run the model
        [S, E, I] = basic_SEI_BestFit(alpha, beta, 1./gamma, 400, [N(j)-k, 0, k]); %#ok<*SAGROW>
        
        % Store the required data
        EIsum = E + I;
        EIdiff = EIsum(:,end) - k;
        EI(k,j) = mean(EIdiff(EIdiff>0));
        
        % Determine the probability of epidemic taking off
        EIprob(k,j) = sum(EIdiff>0)/length(EIdiff);
        
        % Store a matrix of initial conditions
        initialI(k,j) = k;
    end
    
end

% Set nans
initialI(initialI==0)=NaN;

% Plot
figure; set(gcf,'WindowStyle','docked')
subplot(1,2,1)
hand1 = imagesc(EI); set(gca,'YDir','normal'); colorbar
set(hand1,'alphadata',~isnan(initialI))
subplot(1,2,2)
hand2 = imagesc(EIprob); set(gca,'YDir','normal'); colorbar
set(hand2,'alphadata',~isnan(initialI))
