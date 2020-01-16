% Set up the data
%clear; close all
% 
b = 0.01;
alpha = 0.6;
gamma = 0.25; % 1 progression to infectious every 4 days

% alpha = 0.6630;
% b = 0.0047;
% gamma = 3.0608;
tau = 0;

N = 5;

% Create the generator matrix
[Q,HHconfig] = SEI(N); %#ok<*AGROW>

% Set up relevant stuff
beta = b/((N-1)^alpha); %#ok<*PFOUS>
Q.inf = beta*Q.inf;
M = Q.inf + tau*Q.QC + gamma*Q.prog; Mfull = M'; % Transpose the matrix

% Generate the initial conditions vector
P0 = zeros(length(HHconfig.dataI(:,1)),1); P0(length(HHconfig.dataI(:,1))-1,1) = 1;

% Runge-Kutta order 4,5 i.e. ode45
f = @(t,x)Mfull*x;
[time, P] = ode45(f,[0 365],P0);

% Generate the data
for i = 1 : length(P(:,1))
    E(i,1) = (P(i,:)*HHconfig.dataI(:,2));
    S(i,1) = (P(i,:)*HHconfig.dataI(:,1));
    I(i,1) = (P(i,:)*HHconfig.dataI(:,3));
end
hold on
subplot(1,2,1)
set(gcf,'WindowStyle','Docked')
plot(time,S,'k',time,E,'g',time,I,'r'); hand = legend('S','E','I'); set(hand,'Box','off')
xlabel('Time'); ylabel('No of individuals')
box off; hold on
plot(time,S+E+I,'b--'); set(gca,'YLim',[0 11])