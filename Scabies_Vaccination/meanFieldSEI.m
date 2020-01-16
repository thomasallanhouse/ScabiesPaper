% Function contains the system of ODE for the mean field SEI model
% Written by Tim Kinyanjui
% on 8th May 2017

% ODE function
function dy_dt = meanFieldSEI(t0,y0,beta,alpha,gamma) %#ok<INUSL>

% Per capita rate of infection
N = sum(y0);
lambda = beta/((N-1)^alpha);

% ODE functions
dy_dt = zeros(3,1);
% Susceptibles
dy_dt(1) = -lambda*y0(1)*y0(3);

% Latent
dy_dt(2) = lambda*y0(1)*y0(3) - gamma*y0(2);

% Infectious
dy_dt(3) = gamma*y0(2);

return