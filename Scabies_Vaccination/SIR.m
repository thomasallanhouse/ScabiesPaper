function [Qinf, Qrec, Qext, dataI] = SIR(N)

% Purpose:
% Function written to generate the transition rates matrix for household
% model with SIRS disease structure
%
% Input:
%   N: Household size which we hold constant in the simulations
%
% Output:
%   Q: Transition matrix
%   data: Contains the household configurations


% List the household configurations exhaustively i.e. all the possible
% household configurations and exclude the ones that do not retain the
% household size

% Function written by Tim Kinyanjui on 15th Feb 2012
% Inspired by work from J Ross et al.
% Edited on 18th Feb 2012.


% Dimension of the matrix. This only works for 3 epidemiological classes. I
% don't know what will happen when they are not 3, but will find out
% later.
dim = sum(1:(N+1));

Qinf = zeros(dim,dim); %spalloc
Qrec = zeros(dim,dim); %spalloc
Qext = zeros(dim,dim); %spalloc
m=0;
for ss=0:N
    for ii=0:(N-ss)
        m = m+1;
        I{ss+1,ii+1} = m;  %#ok<AGROW>
    end
end
dataI = zeros(dim,3);

% Counter for susceptibles
for ss=0:N
    
    % Counter for infecteds
    for ii=0:(N-ss)
        
        
        % If susceptibles are more than 1 and there is atleast one
        % infected then the rate of infection exists
        if (ss >= 1 && ii >= 1)           
            % Rate of within household infection
            Qinf(I{ss+1,ii+1},I{ss,ii+2}) = ii*ss/(N-1);           
        end
        
        % For recovery to occur, the number of infecteds must be
        % greater than 1
        if ii >= 1            
            % Rate of recovery
            Qrec(I{ss+1,ii+1},I{ss+1,ii}) = ii;           
        end
        if (ss >= 1)           
            % Rate of within household infection
            Qext(I{ss+1,ii+1},I{ss,ii+2}) = ss;           
        end
        
        % Store the relevant indices to help identify the household
        % configurations outside of this function.
        
        dataI(I{ss+1,ii+1},:) = [ss, ii, N-ss-ii];
        
    end
    
end

Qinf=Qinf-diag(sum(Qinf,2));
Qext=Qext-diag(sum(Qext,2));
Qrec=Qrec-diag(sum(Qrec,2));