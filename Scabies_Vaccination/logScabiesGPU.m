function LLhood = logScabiesGPU(params,sMethod,tau,N,I,T)

% Function calculates the likelihood of observing the data given the
% parameters b & alpha for all the carehomes%
% beta (transmission parameter) takes the following form:
%   beta = b/((N-1)^alpha)
%
% Function written by Tim on 6th Aug 2015
% Edited on 16th Sep 2015 t accomodate SEI scabies model
% Univerity of Manchester

% Step size for the Euler backward expansion
h = 1;
Tol = 1e-2;

for k = 1:length(params(:,1))
    
    b = params(k,1);
    alpha = params(k,2);
    gamma = params(k,3);
    
    % If params is negative, reject the calculation and return a big
    % negative infinity
    if b >= 0 && gamma >=0
        
        % Loop through the carehomes data
        parfor i = 1:length(N)
            
            % Create the generator matrix
            [Q,HHconfig] = SEI(N(i)); %#ok<*AGROW>
            
            % Find locations on vector with I(i) infecteds
            locI = (HHconfig.dataI(:,2) + HHconfig.dataI(:,3)) == I(i);
            
            % Set up relevant stuff
            beta(i,1) = b/((N(i)-1)^alpha); %#ok<*PFOUS>
            Q.inf = beta(i,1)*Q.inf;
            M = Q.inf + tau*Q.QC + gamma*Q.prog;
            M = sparse(M');
            
            % Generate the initial conditions vector
            P0 = zeros(length(HHconfig.dataI(:,1)),1); P0(length(HHconfig.dataI(:,1))-1,1) = 1;
            
            % Solve the master equation using the following methods
            if sMethod == 7
                
                % Use RATKRYL
                A = M*T(i);
                [pRat, errest] = sikrylov_proj(A,P0,Tol);
                
                % Normalise
                %pRat = abs(pRat)/sum(pRat);
                
                % Calculate the log-likelihood of observing the data given the parameters
                % llExp(i,1) = pCheb(N(i)+1-I(i));
                locRat = pRat(locI);
                %locPExp = locExp(locExp~=0);
                locPrat = sum(locRat);
                llExp(i,1) = log(locPrat);
                
            elseif sMethod == 1
                
%                 try
%                     
%                 % Solve using Expokit
%                 pExp = mexpv(T(i), M, P0, Tol);
%                 
%                 catch err
%                     
%                     dbstop in logScabiesGPU at 56
%                     
%                 end
                pExp = mexpv(T(i), M, P0, Tol);
                
                % Normalize to guarantee probability vector
                pExp = abs(pExp);
                pExp = pExp/sum(pExp);
                
                % Calculate the log-likelihood of observing the data given the parameters
                % llExp(i,1) = sum(log(pExp(locI)));
                locExp = pExp(locI);
                % locPExp = locExp(locExp~=0);
                locPExp = sum(locExp);
                llExp(i,1) = log(locPExp);
                
            elseif sMethod == 2
                
                % Chebyshev's expansion
                pCheb = polycheby2(M*T(i), P0, Tol, 2850, min(diag(M*T(i))), max(diag(M*T(i))));
                
                % Normalize to guarantee probability vector
                pCheb = abs(pCheb);
                % pCheb(pCheb < 0) = min(pCheb(pCheb > 0));
                pCheb = pCheb/sum(pCheb);
                
                % Calculate the log-likelihood of observing the data given the parameters
                % llExp(i,1) = pCheb(N(i)+1-I(i));
                locExp = pCheb(locI);
                %locPExp = locExp(locExp~=0);
                locPExp = sum(locExp);
                llExp(i,1) = log(locPExp);
                
            elseif sMethod == 3
                
                % Runge-Kutta order 4,5 i.e. ode45
                f = @(t,x)M*x;
                [~, p45_1] = ode45(f,[0 T(i)],P0,odeset('RelTol',Tol));
                p45 = p45_1(end,:)'; % Take the last time point of interest
                
                % Normalize to guarantee probability vector
                p45 = abs(p45);
                p45 = p45/sum(p45);
                
                % Calculate the log-likelihood of observing the data given the parameters
                % llExp(i,1) = p45(N(i)+1-I(i));
                locExp = p45(locI);
                %locPExp = locExp(locExp~=0);
                locPExp = sum(locExp);
                llExp(i,1) = (log(locPExp));
                
                
            elseif sMethod >= 4 && sMethod <= 6
                
                % DA process
                order = sMethod - 3;
                II = eye(length(HHconfig.dataI(:,1)),length(HHconfig.dataI(:,1)));
                pDA_1 = daMethod(h,II,M,T(i),order,P0,Tol);
                pDA = pDA_1(end,:)';
                
                % Normalize to guarantee probability vector
                pDA = abs(pDA);
                pDA = pDA/sum(pDA);
                
                % Calculate the log-likelihood of observing the data given the parameters
                % llExp(i,1) = pDA(N(i)+1-I(i));
                locExp = pDA(locI);
                %locPExp = locExp(locExp~=0);
                locPExp = sum(locExp);
                llExp(i,1) = log(locPExp);
                
            else
                
                error('Please specify a valid method')
                
            end
            
        end
        
        % Calculate the log-likelihood
        LLhood(k,1) = sum(llExp);
        
    else
        
        % Return negative infinity if negative
        LLhood(k,1) = -inf;
        warning('Input parameter negative')
        
    end
    
end

return