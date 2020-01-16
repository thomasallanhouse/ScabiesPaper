function [LLhood, timeD, aGrid, bGrid, gGrid] = runScabies(startMethod, endMethod)

% Function runs the Scabies model
% Input takes a scalar that corresponds to the following:
% 1: Expokit
% 2: Chebyshev expansion
% 3: Runge-Kutta (4,5) i.e. ode45
% 4: DA order 1
% 5: DA order 2
% 6: DA order 3 ... an so on
%
% Written by Tim on 16th Sept 2015 to accomodate SEI scabie model
% Univesity of Manchester

% Load the data
data = importdata('nit.txt');

% Get data in the proper vectors
[N, I, T] = deal(data(:,1), data(:,2), data(:,3));

% Set-up grid
aGrid = linspace(1e-3,1.2,2);
bGrid = linspace(1e-3,0.008,2);
gGrid = linspace(1,6,2);


% aGrid = 0:1:2;
% bGrid = 0:0.025:0.05;

% Optimal using fmincon
% aGrid = 0.6630;
% bGrid = 0.0047;
% gGrid = 3.0608;

tau = 0;

% Initialise
LLhood = zeros(length(aGrid),length(bGrid));
timeD = zeros(endMethod,1);

% Loop through the different methods
for k = startMethod:endMethod
    
    tic;
    
    for i = 1: length(aGrid)
        
        aaGrid = aGrid(i);
        
        for m = 1 : length(gGrid)
            
            ggGrid = gGrid(m);
            
            parfor j = 1 : length(bGrid)
                
                LLhood(i,j,m,k) = logScabiesGPU([bGrid(j),aaGrid,ggGrid],k,tau,N,I,T);
                
            end
            
        end
        
    end
    
    % Total time it takes to populate the uniform grid
    timeD(k,1) = toc;
    
end

% figure; set(gcf,'WindowStyle','Docked')
% for k = startMethod:endMethod
%     subtightplot(2,3,k,[0.09,0.04],[0.09 0.05],[0.05 0.05])
%     messT = sprintf('Time: %.3f',timeD(k,1));
%     contour(bGrid,aGrid,exp(real(LLhood(:,:,k))))
%     title(messT)
%     colormap jet; box on; % colorbar
%     if k==1 || k==4
%         ylabel('Damping with size, \alpha')
%     end
%     if k>=4
%         xlabel({'Transmissibility, \beta (days^{-1})'})
%     end
% end

return