% Load the data
data = importdata('nit.txt');

% Get data in the proper vectors
[N, I, T] = deal(data(:,1), data(:,2), data(:,3));

% Anonymous function to calculate the negative log-likelihood
llhood = @(params)-logScabiesGPU([0.0047,0.663,params],3,0,N,I,T);

% Initial guess
p0 = 3;

% Set the optimization options - use parallel
options = optimset('fmincon');
options.UseParallel = 1;
options.TolX = 1e-10;

% Run the optimisation
tic
pSol = fmincon(llhood,p0,[],[],[],[],0,[],[],options);
time = toc;

% p = linspace(0,20,20);
% parfor i=1:20
%     LLhood(i,1) = llhood(p(i));
% end