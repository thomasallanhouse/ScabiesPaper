% Prepare the workspace
clearvars; %close all

% Load the data
% load multipleChains_041016.mat
% load multipleChains_081217
% load multipleChains_040417
% load multipleChains_200417
% load multipleChains_240417
% load multipleChains_020517
% load multipleChains_150517
% load multipleChains_190517
% load multipleChains_200917
% load multipleChains_260917
load multipleChains_021017
% load multipleChains_051017

% Calculate the Effective Sample Size
burnin = 10000;
thin = 1;
chainBB = chainB(burnin:thin:end,1,:); chainBB = chainBB(:);
chainA = chainB(burnin:thin:end,2,:); chainA = chainA(:);
chainG = chainB(burnin:thin:end,3,:); chainG = chainG(:);
chainLL = chain_ll(burnin:thin:end,:); chainLL = chainLL(:);

[acf,lags] = autocorr(chainBB,20);
ESS_beta = length(chainBB)/(1 + (2*(sum(acf))));

[acf,lags] = autocorr(chainA,20);
ESS_alpha = length(chainA)/(1 + (2*(sum(acf))));

[acf,lags] = autocorr(chainG,20);
ESS_gamma = length(chainG)/(1 + (2*(sum(acf))));




% Burn-in and thin
thin = 20;
boundA = 800;
boundB = 4000;
boundG = 300;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot the autocorrelation function
figure; set(gcf,'WindowStyle','docked')
subplot(1,3,1)
BB = chainB(burnin:end,1,:); BB = BB(:);
[acf,lags,bounds] = autocorr(BB,boundA);
plot(lags,acf,'r'); hold on
plot([min(lags) max(lags)],[bounds(1) bounds(1)],'b','LineStyle','--')
plot([min(lags) max(lags)],[bounds(2) bounds(2)],'b','LineStyle','--')
plot([min(lags) max(lags)],[0 0],'k','LineStyle','-')
title('Auto correlation Function - \beta')
xlabel('Lag'); ylabel('Sample Autocorrelation')
box off

subplot(1,3,2)
AA = chainB(burnin:end,2,:); AA = AA(:);
[acf,lags,bounds] = autocorr(AA,boundB);
plot(lags,acf,'r'); hold on
plot([min(lags) max(lags)],[bounds(1) bounds(1)],'b','LineStyle','--')
plot([min(lags) max(lags)],[bounds(2) bounds(2)],'b','LineStyle','--')
plot([min(lags) max(lags)],[0 0],'k','LineStyle','-')
%set(gca,'XLim',[0 255])
title('Auto correlation Function - \alpha')
xlabel('Lag'); ylabel('Sample Autocorrelation')
box off

subplot(1,3,3)
GG = chainB(burnin:end,3,:); GG = GG(:);
[acf,lags,bounds] = autocorr(GG,boundG);
plot(lags,acf,'r'); hold on
plot([min(lags) max(lags)],[bounds(1) bounds(1)],'b','LineStyle','--')
plot([min(lags) max(lags)],[bounds(2) bounds(2)],'b','LineStyle','--')
plot([min(lags) max(lags)],[0 0],'k','LineStyle','-')
title('Auto correlation Function - \gamma')
xlabel('Lag'); ylabel('Sample Autocorrelation')
box off

% Burn-in and thin
chainBB = chainB(burnin:thin:end,1,:); chainBB = chainBB(:);
chainA = chainB(burnin:thin:end,2,:); chainA = chainA(:);
chainG = chainB(burnin:thin:end,3,:); chainG = chainG(:);
chainLL = chain_ll(burnin:thin:end,:); chainLL = chainLL(:);

% Marginals
figure; set(gcf,'WindowStyle','Docked')

subtightplot(3,3,1,[0.07 0.03],[0.1 0.05],0.055)
hand = histogram(chainBB,30,'Normalization','pdf'); hold on
set(hand,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5])
xx = 0:0.0001:0.06;
plot(xx,exppdf(xx,1/10),'k','LineWidth',2)
% muRand = exprnd(1/10^2,1000000,1);
% [counts,cent] = hist(muRand,20); counts = counts/sum(counts);
% plot(cent,counts,'k','LineWidth',2)
box off; set(gca,'XLim',[0 0.05]) % 0.04

subtightplot(3,3,4,[0.07 0.03],[0.1 0.05],0.055)
% scatter(chainBB,chainA,4,(chainLL),'o','filled')
%[N,C] = hist3([chainBB(:), chainA(:)],[80 80]);
%contour(C{1},C{2},N',15);
[~,density,X,Y]=kde2d([chainBB chainA],2000,[0 -1],[0.06 1.5]);
contour(X,Y,density,20)
box off; yhand = ylabel('\alpha'); set(yhand,'FontSize',14)
set(gca,'XLim',[0 0.03])

subtightplot(3,3,5,[0.07 0.03],[0.1 0.05],0.055)
hand = histogram(chainA,20,'Normalization','pdf'); hold on
set(hand,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5])
plot([-1 1.5],0.5*[1 1],'k','LineWidth',2)
% muA = -1 + (1+1)*rand(1000000,1);
% [counts,cent] = hist(muA,20); counts = counts/sum(counts);
% plot(cent,counts,'k','LineWidth',2)
box off; %colormap jet; 
set(gca,'XLim',[-1 1.5])

subtightplot(3,3,7,[0.07 0.03],[0.115 0.05],0.05)
% scatter(chainBB,chainG,4,(chainLL),'o','filled')
%[N,C] = hist3([chainBB(:), chainG(:)],[50 50]);
%contour(C{1},C{2},N',15);
[~,density,X,Y]=kde2d([chainBB chainG],2000,[0 0.1],[0.028 20]);
contour(X,Y,density,20)
box off; xhand = xlabel('\beta'); yhand = ylabel('\gamma'); %colormap jet;
set([xhand yhand],'FontSize',14); set(gca,'XLim',[0 0.03])

subtightplot(3,3,8,[0.07 0.03],[0.115 0.05],0.05)
% scatter(chainA,chainG,4,(chainLL),'o','filled')
%[N,C] = hist3([chainA(:), chainG(:)],[50 50]);
%contour(C{1},C{2},N',15);
[~,density,X,Y]=kde2d([chainA chainG],2000,[-1 0.1],[1.5 20]);
contour(X,Y,density,20)
box off; xhand = xlabel('\alpha'); set(xhand,'FontSize',14)
set(gca,'XLim',[-1 1.5])

subtightplot(3,3,9,[0.07 0.03],[0.115 0.05],0.05)
%hist(chainG,16)
hand = histogram(chainG,20,'Normalization','pdf'); hold on
set(hand,'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5])
data1 = importdata('gammaFitdata.mat');
xx = 0:0.01:30;
plot(xx,gampdf(xx,data1.xopt(1),data1.xopt(2)),'k','LineWidth',2)
% muG = gamrnd(data1.xopt(1),data1.xopt(2),[1000000,1]);
% [counts,cent] = hist(muG,20); counts = counts/sum(counts);
% plot(cent,counts,'k','LineWidth',2)
box off; xhand = xlabel('\gamma'); set(xhand,'FontSize',14)
colormap bone

% % 3D scatter plot
% figure; set(gcf,'WindowStyle','Docked')
% scatter3(chainBB,chainA,chainG,6,exp(chainLL),'o','filled')
% colorbar; %colormap jet
% xlabel('\beta'); ylabel('\alpha'); zlabel('\gamma')

% Figure to show convergence and mixing of chains
figure; set(gcf,'WindowStyle','docked')
subtightplot(3,1,1,[0.05 0.1])
plot(1:length(chainBB),chainBB)
hold on
%plot([0 length(chainTime)],[0.0047 0.0047],'r')
ylabel('\beta'); set(gca,'XTickLabel',''); box off

subtightplot(3,1,2,[0.05 0.1])
plot(1:length(chainBB),chainA)
hold on
%plot([0 length(chainTime)],[0.663 0.663],'r')
ylabel('\alpha'); set(gca,'XTickLabel','','YTick',0:0.25:1); box off

subtightplot(3,1,3,[0.05 0.1])
plot(1:length(chainBB),chainG)
hold on
%plot([0 length(chainTime)],[1/3.0608 1/3.0608],'r')
ylabel('\gamma'); xlabel('Chain length'); box off

