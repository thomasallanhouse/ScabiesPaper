function [Q, HHconfig] = SEIR(N)

% Function written to generate the transition rates matrix for household
% model with SEIR disease structure
%
% Input:
%   N: Household size which we hold constant in the simulations
%
% Output:
%   Q: Transition matrix
%   dataI: Contains the household configurations
%
% Function written by Tim Kinyanjui on 16th Sept 2015
% University of Manchester


% List the household configurations exhaustively i.e. all the possible
% household configurations and exclude the ones that do not retain the
% household size

m = 0;

% Transforms two integers to locations where we will store the variables
for ss=0:N
    
    for ii=0:(N-ss)
        
        for jj = 0:(N-ss-ii)
            
            m = m+1;
            
            I{ss+1,ii+1,jj+1} = m;  %#ok<AGROW>
            
        end
        
    end
    
end

% Determine the dimension of the system. Otherwise, the matrices will not
% be square

dim = m;

Q.inf = spalloc(dim,dim,3*dim);
Q.QC = Q.inf;
Q.prog = Q.inf;
Q.Vacc = Q.inf;

dataI = zeros(dim,4);

% Counter for susceptibles
for ss=0:N
    
    for ii=0:(N-ss)
        
        for jj = 0:(N-ss-ii)
            
            % If susceptibles are less than N then there are infecteds and therefore
            % infection within the household can happen
            if  ss > 0 && jj > 0
                
                % Rate of within household infection
                % Q.inf(I{ss+1,ii+1},I{ss}) = betaHH*ss*(N-ss)/(N-1);
                Q.inf(I{ss+1,ii+1,jj+1},I{ss+1-1,ii+1+1,jj+1}) = ss*jj;
                
            end
            
            % Infection from outside of the household
            % Infection will only occur if number of susceptibles is greater than 0
            if ss > 0
                
                % Infection from outside of the household
                % Q.QC(I{ss+1},I{ss+1-1}) = ss;
                Q.QC(I{ss+1,ii+1,jj+1},I{ss+1-1,ii+1+1,jj+1}) = ss;
                
            end
            
            % For progression to active infection to occur
            if ii >= 1
                
                % Progression to infection
                Q.prog(I{ss+1,ii+1,jj+1},I{ss+1,ii+1-1,jj+1+1}) = ii;
                
            end
            
            % Vaccination
            if jj >= 1
                
                % Vaccinaton happens
                Q.Vacc(I{ss+1,ii+1,jj+1},I{ss+1,ii+1,jj+1-1}) = jj;
                
            end
            
            % Store the relevant indices to help identify the household
            % configurations outside of this function.
            
            dataI(I{ss+1,ii+1,jj+1},:) = [ss, ii, jj, N-ss-ii-jj];
            
        end
        
    end
    
end

Q.inf = Q.inf-diag(sum(Q.inf,2));
Q.QC = Q.QC-diag(sum(Q.QC,2));
Q.prog = Q.prog-diag(sum(Q.prog,2));
Q.Vacc = Q.Vacc-diag(sum(Q.Vacc,2));

HHconfig.dataI = dataI;

return