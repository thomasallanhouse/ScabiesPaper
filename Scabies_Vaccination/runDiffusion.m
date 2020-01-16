% Script runs the stochastic SEI model
% Written by Tim Kinyanjui
% 16th June 2016

% Clear workspace
clearvars;

% Determine the parameters to use from the MCMC results
% load multipleChains_081217
load multipleChains_021017.mat
burnin = 10000;
thin = 20;
chainBB = chainB(burnin:thin:end,1,:); beta = chainBB(:);
chainA = chainB(burnin:thin:end,2,:); alpha = chainA(:);
chainG = chainB(burnin:thin:end,3,:); gamma = chainG(:);

% Run the diffusion model
parfor i = 1:30
    [S(:,:,i), E(:,:,i), I(:,:,i)] = basic_SEI(alpha, beta, 1./gamma);
end

% Take the mean of the replicates
NN = S+E+I; N = NN(1);
SS = mean(S,3); EE = mean(E,3); II = mean(I,3); EI = (E+I)/N; EEII = mean(EI,3);


% Do the plotting
% S
figure
set(gcf,'WindowStyle','docked')
%plot(0:365,SS,'Color',[0.5 0.5 0.5])
hold on
plot([0:365]./7,mean(SS,1),'k','LineWidth',2)
plot([0:365]./7,prctile(SS,[2.5 97.5],1),'k--','LineWidth',2)
title('S'); set(gca,'XTick',0:2:52,'XLim',[0 52])

% E
figure
set(gcf,'WindowStyle','docked')
%plot(0:365,EE,'b')
hold on
plot([0:365]./7,mean(EE,1),'k','LineWidth',2)
plot([0:365]./7,prctile(EE,[2.5 97.5],1),'k--','LineWidth',2)
title('E'); set(gca,'XTick',0:2:52,'XLim',[0 52])

% I
figure
set(gcf,'WindowStyle','docked')
%plot(0:365,II,'Color',[0.5 0.5 0.5])
hold on
plot([0:365]./7,mean(II,1),'k','LineWidth',2)
plot([0:365]./7,prctile(II,[2.5 97.5],1),'k--','LineWidth',2)
title('I'); set(gca,'XTick',0:2:52,'XLim',[0 52])

% Calculate the area under the curve using Trapezoidal rule
Iq = II'/N;
for i = 1:length(Iq(1,:))
    for j = 1:(length(Iq(:,1))-1)
        vecA = Iq(j:j+1,i);
        auc(j,i) = trapz(vecA);
    end
end
cumAuc = cumsum(auc,1); % Take along the columns

figure;
set(gcf,'WindowStyle','docked')
subplot(1,2,1)
%plot(0:365,auc,'Color',[0.5 0.5 0.5])
hold on
plot([0:365]./7,mean(EEII,1),'k','LineWidth',2)
plot([0:365]./7,prctile(EEII,[2.5 97.5],1),'k--','LineWidth',2)
title('Economic cost'); set(gca,'XTick',0:3:52,'XLim',[0 52]);xlabel('Time in weeks'); ylabel('Proportion of people to treat (E+I)')
hand = text(3,0.9,'A');
set(hand,'FontWeight','bold')

subplot(1,2,2)
hold on
plot([0:364]./7,mean(cumAuc,2),'k','LineWidth',2)
plot([0:364]./7,prctile(cumAuc,[2.5 97.5],2),'k--','LineWidth',2)
title('QALY cost'); set(gca,'XTick',0:3:52,'XLim',[0 52]);xlabel('Time in weeks'); ylabel('QALY cost')
hand = text(3,265,'B');
set(hand,'FontWeight','bold')