% Function runs the diffusion (Langevin) approximation to the SEI model
% alpha, beta  and gamma can be vectors and the code is vectorised
% Written by Tim Kinyanjui on 16th June 2016
function [S, E, I] = basic_SEI(alpha, beta, gamma)

% Change input to be vectors of length n x 1
beta = beta(:);
gamma = gamma(:);
alpha = alpha(:);

% Initialize the parameters and preallocate
t = 365;
S = zeros(size(beta, 1), t);
E = S;
I = S;

% Initial state conditions
%mu = .05/365;
mu = 0;
S(:, 1) = 45;
E(:, 1) = 0;
I(:, 1) = 5;
N = S(:,1)+E(:,1)+I(:,1);

for curr_t = 1:t
    
    % Determine the demogrpahic events. Note the diffusion approx
    newbirths = round(max(0, mu.*N + sqrt(mu.*N)*randn(1, 1)));
    deaths = round(max(0, mu.*N + sqrt(mu.*N)*randn(1, 1)));
    Sdeaths = floor(S(:,curr_t)./N.*deaths);
    Edeaths = floor(E(:,curr_t)./N.*deaths);
    Ideaths = floor(I(:,curr_t)./N.*deaths);
    
    % Total population
    N = N + newbirths - Sdeaths - Edeaths - Ideaths;
    
    % Determine the epidemiological events
    % Infection
    mu_inf = S(:, curr_t).*beta.*I(:,curr_t)./((N-1).^alpha);
    newinf = round(min(S(:, curr_t), max(0, mu_inf + sqrt(mu_inf).*randn(size(mu_inf, 1), 1))));
    
    % Progression to infected from E state
    mu_rec =gamma.*E(:,curr_t);
    newrec = round(min(E(:, curr_t), max(0, mu_rec + sqrt(mu_rec).*randn(size(mu_rec, 1), 1))));
    
    % Update the system depending on the events
    S(:, curr_t+1) = S(:, curr_t) + newbirths - newinf - Sdeaths;
    E(:, curr_t+1) = E(:, curr_t) + newinf - newrec - Edeaths;
    I(:, curr_t+1) = I(:, curr_t) + newrec - Ideaths;
    
end

% End of function
return