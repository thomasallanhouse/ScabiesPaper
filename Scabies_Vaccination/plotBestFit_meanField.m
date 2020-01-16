
% Plot the best fit
% Written by Tim Kinyanjui
% on 8th May 2017
% University of Manchester
 
% Clear workspace
clearvars;
timmer = tic;
 
% Determine the parameters to use from the MCMC results
%load multipleChains_081217
% load multipleChains_150517
% load multipleChains_190517.mat
load multipleChains_021017.mat

burnin = 10000;
thin = 20;
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
 
% Run the mean field model
 
for j = 1:length(N)
   
    %[S(j,kk).SS, E(j,kk).EE, I(j,kk).II] = runMeanFieldSEI(alpha, beta, 1./gamma, T(j), [N(j)-1, 0, 1]); %#ok<SAGROW>
    % [time, Y(:,:,:,j)] = runMeanFieldSEI(linspace(0,T(j),100),[N(j)-1, 0, 1],beta,alpha,1./gamma);
    [Y(:,j), fullSet(:,j)] = SEI_Gillespie_multiple(370,[N(j)-1, 0, 1],beta,alpha,gamma);
   
end
 
% Do the plotting
hand1 = figure; set(hand1,'WindowStyle','Docked')
 
k=1;
for j = 1:length(N)
   
    %%%%%%%%%%%%%%%%%%%%%%%%% Plot the model fit for rach household
    figure(hand1)
    subtightplot(7,2,k,[0.01 0.05],[0.1 0.05],0.055)
    for i = 1:length(fullSet(:,1))
        % plot(fullSet{i,j}(:,1),sum(fullSet{i,j}(:,3:end),2),'Color',[0.5 0.5 0.5]); hold on
        Is(:,i) = sum(fullSet{i,j}(:,3:end),2);
        justIs(:,i) = fullSet{i,j}(:,4);
    end
    
    timeIs = fullSet{1,j}(:,1);
    plot(timeIs,mean(Is,2),'k','LineWidth',2); hold on
    plot(timeIs,prctile(Is,[2.5 97.5],2),'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','--')
    scatter(T(j),II(j),'bo','filled')
    xxlim = get(gca,'Xlim');
    plot(xxlim,[N(j) N(j)],'r','LineWidth',2); box off;
    if j ~= length(N)
        set(gca,'XTickLabel','')
    end
    axis tight
   
    %%%%%%%%%%%%%%%%%%%%%%%%% Plot the QALY each household
   
    % Calculate the area under the curve using Trapezoidal rule
    Iq = justIs/N(j);
    for i = 1:length(Iq(1,:))
        for kk = 1:(length(Iq(:,1))-1)
            vecA = Iq(kk:kk+1,i);
            auc(kk,i) = trapz(vecA);
        end
    end
    
    cumAuc = cumsum(auc,1); % Take along the columns
    subtightplot(7,2,k+1,[0.01 0.05],[0.1 0.05],0.055)
    hold on
    plot(timeIs(1:end-1),mean(cumAuc,2),'k','LineWidth',2)
    plot(timeIs(1:end-1),prctile(cumAuc,[2.5 97.5],2),'LineWidth',2,'Color',[0.5 0.5 0.5],'LineStyle','--')
    if j ~= length(N)
        set(gca,'XTickLabel','')
    end
    axis tight
    %title('QALY cost'); set(gca,'XTick',0:3:52,'XLim',[0 52]);xlabel('Time in weeks'); ylabel('QALY cost')
    %hand = text(3,265,'B');
    %set(hand,'FontWeight','bold')
   
%     %%%%%%%%%%%%%%%%%%%%%%%%% Plot the Economic cost for each household
%     subtightplot(7,3,k+1,[0.01 0.05],[0.1 0.05],0.055)
%     plot(timeIs,mean(Is,2)/N(j),'k','LineWidth',1); hold on
%     plot(timeIs,prctile(Is,[2.5 97.5],2)/N(j),'k','LineWidth',1,'Color',[0.5 0.5 0.5])
%     % title('Economic cost'); set(gca,'XTick',0:3:52,'XLim',[0 52]);xlabel('Time in weeks'); ylabel('Proportion of people to treat (E+I)')
%     % hand = text(3,0.9,'A');
%     if j ~= length(N)
%         set(gca,'XTickLabel','')
%     end
%     axis tight; box off
   
    % Increase the subplot counter
    k = k + 2;
 
end

% Store the computational time
timer = toc(timmer);

% Save the results
save Figure3Results_121017