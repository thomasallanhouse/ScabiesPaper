
function [time1, Y] = SEIR_Gillespie(sim_time,Y0,beta,alpha,gamma,tau)
 
% This function implements the stochastic SEI model without demography.
%  two independent events can evolve the system
%   1) Infection from S->E
%   2) Progression from E to I E->I
%   3) Vaccination from I to R
% Implemented using the Gillespie's First reaction algorithm.
% Input:
%   sim_time: Vector with thSimulation time
%   Y0: Iniial conditions [S, E, I]
%   beta: Transmission parameter
%   alpha: Contact scaling
%   gamma: Rate of progression to infectious
%
% Function written and designed by Tim kinyanjui
% 8th May 2017
% University of Manchester - School of Mathematics
% % Initialize the vectors for code optimization
% S=zeros(sim_time,1);
% I=zeros(sim_time,1);
% time=zeros(sim_time,1);
% rate=zeros(5,1);
% dt=zeros(5,1);
% %  End of initialization
% Store the initial conditions

S(1) = Y0(1);
E(1) = Y0(2);
I(1) = Y0(3);
R(1) = Y0(4);
time(1) = 0;

% Set up data
N = sum(Y0);
lambda = beta/((N-1)^alpha);

% Initilaize counter to keep track of vector positions
i = 2;

% Now implement the G-Algorithm
while(time(i-1) <= sim_time)
    
    if  E(i-1) > 0 || (S(i-1) > 0 && I(i-1) > 0) || I(i-1) > 0
        
        % Determine the Rates at which events will occur
        
        % Rate of infection
        rate(1,1) = lambda*S(i-1)*I(i-1);
        
        % Rate of progression to infectious state
        rate(2,1) = gamma*E(i-1);
        
        % Rate of vaccination
        rate(3,1) = tau*I(i-1);
        
        % Determine the times to the next event for each event (dt)
        dt = -log(rand(3,1))./rate;
        
        % Find the minimum of the dt to determine the time to the next
        % event and the next event is determined by the position (on the
        % vector) of the minimum value returned by the variable k
        [~, k] = min(dt);
        
        % Update the system states
        if(k == 1)
            
            % Infection occurs
            S(i,1) = S(i-1) - 1;
            E(i,1) = E(i-1) + 1;
            I(i,1) = I(i-1);
            R(i,1) = R(i-1);
            
        elseif(k == 2)
            
            % Progression to infectious
            S(i,1) = S(i-1);
            E(i,1) = E(i-1) - 1;
            I(i,1) = I(i-1) + 1;
            R(i,1) = R(i-1);
            
        elseif(k==3)
            
            % Vaccination
            S(i,1) = S(i-1);
            E(i,1) = E(i-1);
            I(i,1) = I(i-1) - 1;
            R(i,1) = R(i-1) + 1;
            
        else
            
            warning('Nothing happening: Should really not be here')
            
        end
        
        % Progress time by dt(k)
        time(i) = time(i-1) + dt(k);
        i=i+1;
        
        % For debugging
        if sum(isinf(time)) > 0
            dbstop in SEIR_Gillespie at 114
        end
        
        %elseif S(i-1) == 0 && E(i-1) == 0 && I(i-1) == 0
    else
        
        % Stop and the state of the system remains the same
        S(i,1) = S(i-1);
        E(i,1) = E(i-1);
        I(i,1) = I(i-1);
        R(i,1) = R(i-1);
        time(i) = sim_time;
        
        % Get out of the while loop (Not good programming practice)
        break;
        
    end
    
    %    end
    
end

% Final data
Y1 = [S E I R];

% Interpolate the output at the same time points
time1 = (0:1:sim_time)';
% if sum(isinf(time)) > 0
%     dbstop in SEI_Gillespie at 120
% end
Y(:,1) = interp1(time,Y1(:,1),time1);
Y(:,2) = interp1(time,Y1(:,2),time1);
Y(:,3) = interp1(time,Y1(:,3),time1);
Y(:,4) = interp1(time,Y1(:,4),time1);

return