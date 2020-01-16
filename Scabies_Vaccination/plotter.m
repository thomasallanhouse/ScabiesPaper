% Script plots the results of the mcmc runs from the scabies model
% Written by Tim Kinyanjui on the 17th of Oct 2015
% University of Manchester

% Prepare the workspace
clearvars

% Load the data
load results090616.mat


%load unifPriorsTemp

% Burn-in and thin
thin = 10;
burnin = 3000;
chainBB = chainB(burnin:thin:end,1);
chainA = chainB(burnin:thin:end,2);
chainG = chainB(burnin:thin:end,3);
chainLL = chain_ll(burnin:thin:end,1);
chainTimeThind = chainTime(burnin:thin:end);

% Figure to show convergence and mixing of chains
figure; set(gcf,'WindowStyle','docked')
subtightplot(3,1,1,[0.05 0.1])
plot(chainTime,chainB(:,1))
hold on
%plot([0 length(chainTime)],[0.0047 0.0047],'r')
ylabel('\beta'); set(gca,'XTickLabel',''); box off

subtightplot(3,1,2,[0.05 0.1])
plot(chainTime,chainB(:,2))
hold on
%plot([0 length(chainTime)],[0.663 0.663],'r')
ylabel('\alpha'); set(gca,'XTickLabel','','YTick',0:0.25:1); box off

subtightplot(3,1,3,[0.05 0.1])
plot(chainTime,chainB(:,3))
hold on
%plot([0 length(chainTime)],[1/3.0608 1/3.0608],'r')
ylabel('\gamma'); xlabel('Chain length'); box off

% Marginals
figure; set(gcf,'WindowStyle','Docked')

subtightplot(3,3,1,[0.05 0.03],[0.1 0.05],0.055)
hist(chainBB,6)
box off

subtightplot(3,3,4,[0.05 0.03],[0.1 0.05],0.055)
scatter(chainBB,chainA,4,exp(chainLL),'o','filled')
box off; ylabel('\alpha'); %colormap jet

subtightplot(3,3,5,[0.05 0.03],[0.1 0.05],0.055)
hist(chainA,6)
box off; %colormap jet; 

subtightplot(3,3,7,[0.05 0.03],[0.1 0.05],0.05)
scatter(chainBB,chainG,4,exp(chainLL),'o','filled')
box off; xlabel('\beta'); ylabel('\gamma'); %colormap jet;

subtightplot(3,3,8,[0.05 0.03],[0.1 0.05],0.05)
scatter(chainA,chainG,4,exp(chainLL),'o','filled')
box off; xlabel('\alpha'); %colormap jet

subtightplot(3,3,9,[0.05 0.03],[0.1 0.05],0.05)
hist(chainG,6)
box off; xlabel('\gamma'); %colormap jet

% 3D scatter plot
figure; set(gcf,'WindowStyle','Docked')
scatter3(chainBB,chainA,chainG,6,(chainLL),'o','filled')
colorbar; %colormap jet
xlabel('\beta'); ylabel('\alpha'); zlabel('\gamma')