% Plot the best fit
% Written by Tim Kinyanjui
% on 26th April 2017
% University of Manchester

% Clear workspace
clearvars;

% Determine the parameters to use from the MCMC results
% load multipleChains_081217
load multipleChains_240417
burnin = 1000;
thin = 15;
chainBB = chainB(burnin:thin:end,1,:); beta = chainBB(:);
chainA = chainB(burnin:thin:end,2,:); alpha = chainA(:);
chainG = chainB(burnin:thin:end,3,:); gamma = chainG(:);

% Load the data from the care homes
data = importdata('nit.txt');
[N, II, T] = deal(data(:,1), data(:,2), data(:,3));

% Run the diffusion model
% parfor i = 1:30
%     [S(:,:,i), E(:,:,i), I(:,:,i)] = basic_SEI(alpha, beta, 1./gamma);
% end

% Run the Gillespie model
for kk = 1:1
    
    for j = 1:length(N)
        
        [S(j,kk).SS, E(j,kk).EE, I(j,kk).II] = basic_SEI_BestFit(alpha, beta, 1./gamma, T(j), [N(j)-1, 0, 1]); %#ok<SAGROW>
        
    end

end

% Do the plotting
figure; set(gcf,'WindowStyle','Docked')

for j = 1:length(N)
    
    for kk = 1:1
        % Extract the final Is from the runs
        Iend1(:,kk) = I(j,kk).II(:,end);
    end
    
    % Mean of the repeats
    Iend = mean(Iend1,2);
    
    % Plot
    subtightplot(3,3,j,[0.07 0.03],[0.1 0.05],0.055)
    hand = histogram(Iend,15); 
    % [counts, cent] = hist(Iend,15); plot(cent, counts,'r');
    hold on; yylim = get(gca,'YLim');
    set(hand,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5])
    plot([II(j) II(j)],[yylim(1) yylim(2)],'LineWidth',2,'Color','k'); box off
    set(gca,'YLim',yylim)

end
