%Script to run other algorithms
addpath(genpath('../AMP'))
addpath('../algorithm')
%addpath('../outdated')

optIn.tryBigampEM = 0;
optIn.tryKsvd = 0;
optIn.tryErspud = 0;
optIn.trySpams = 1;
optIn.tryBCD1 = 1;
optIn.tryBCD2 = 1;
%Specify maximum allowed trials
optIn.maxTrials = 3;
%Problem dimensions
optIn.M = 10; %size of signal
optIn.N = optIn.M; %size of dictionary
optIn.pos = 0;% set coef to be positive
optIn.L = ceil(40*optIn.N); %From Spielman et al.

%SNR
optIn.SNR = 100;

%Specify coding mode (0 for OMP, 1 for TST)
%The specified function is used to compute a coding error for the
%resulting dictionary for each algorithm. TST may be significantly
%slower.
optIn.useTST = 0;

%Precondition option (enable to use a preconditioner for EM-BiG-AMP)
optIn.precondition = 0;

M_list = 10:2:20;
K_list = 1:1:10;
y = cell(length(K_list), length(M_list));
rng(mod(datenum(now)*2^32,2^32))% set seed to be different across time
ttt = tic;
for i = 1:length(K_list)
    for j = 1:length(M_list)
        optIn.K = K_list(i);%Specify sparsity
        optIn.M = M_list(j);
        optIn.N = optIn.M;
        y{i,j} = trial_DL_new(optIn);
    end
end
duration = toc(ttt);
save(strcat('sim4_',num2str(randi(1e8)),'.mat'),'y','K_list','mu_list','optIn', 'duration')
