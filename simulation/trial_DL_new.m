function results = trial_DL_new(optIn)

%trial_DL: This function runs several algorithms on a sample
%instance of the dictionary learning problem. The function can be run with no
%arguments. See the nargin==0 block below for the required fields to call
%the function with optional parameters. The function returns a structure of
%algorithm results.
%This code can be used to produce the results for the noise free phase
%plane plots for dictionary learning in the BiG-AMP arXiv submission. 
%(This code conducts a single trial at 1 point in the phase plane.)

%Add paths
setup_DL

%% Problem Setup

%Turn algs on/off
tryBigampEM = optIn.tryBigampEM;
tryKsvd = optIn.tryKsvd;
tryErspud = optIn.tryErspud;
trySpams = optIn.trySpams;
tryBCD1 = optIn.tryBCD1;
tryBCD2 = optIn.tryBCD2;

%SNR
SNR = optIn.SNR;

%Coding
useTST = optIn.useTST;

%Precondition
precondition = optIn.precondition;

%Max trials
maxTrials = optIn.maxTrials;

%Define problem dimensions
M = optIn.M; % dimension
L = optIn.L; % sample size
N = optIn.N; % equal to M
K = optIn.K; % number of non-zero entries in coef

%Set options
opt = BiGAMPOpt; %initialize the options object
opt.nit = 500; %limit iterations

%Set sizes
problem = BiGAMPProblem();
problem.M = M;
problem.N = N;
problem.L = L;


%% Build the dictionary

%Draw randomly
if false
    mu = optIn.mu;
    A = (ones(M,N)*mu + (1 - mu)*eye(M))^.5; %WARN: I have changed this part.
else
    A = randn(M, N); % default is to use random ref dictionary
end

%Normalize the columns
A = A*diag(1 ./ sqrt(abs(diag(A'*A))));

%Dictionary error function
dictionary_error_function =...
    @(q) 20*log10(norm(A -...
    q*find_permutation(A,q),'fro')/norm(A,'fro'));



%% Compute coefficient vectors

%Compute true vectors with exactly K non-zeroes in each
if optIn.pos == 1
    X = abs(randn(N,L)); %WARN: I changed this to be non-negative !
else
    X = randn(N,L); % default
end
for ll = 1:L
    yada = randperm(N);
    yada2 = zeros(N,1);
    yada2(yada(1:K)) = 1;
    X(:,ll) = X(:,ll) .* yada2; % coef is sparse Gaussian.
end



%% Form the output channel

%Compute noise free output
Z = A*X; % Z is the noise free signal.

%Define the error function
error_function = @(qval) 20*log10(norm(qval - Z,'fro') / norm(Z,'fro'));
opt.error_function = error_function;


%Determine nuw
nuw = norm(reshape(Z,[],1))^2/M/L*10^(-SNR/10);

%Noisy output channel
Y = Z + sqrt(nuw)*randn(size(Z));

%Coding error
coding_error_function = @(q) 20*log10(coding_error(Y,q,K,useTST));


%Initialize results as empty
results = [];


%% EM BiG AMP

if tryBigampEM
    
    %Silence
    opt.verbose = false;
    disp('Starting EM-BiG-AMP')
    
    %Coding
    bestError = inf;
    bestSparsity = inf;
    
    %Compute preconditioner
    if precondition
        Q = chol(inv(Y*Y'));
    else
        Q = 1;
    end
    
    %Define the error function
    QZ = Q*Y; %Using (possibly) noisy data for error function
    error_function2 = @(qval) 20*log10(norm(qval - QZ,'fro') / norm(QZ,'fro'));
    opt.error_function = error_function2;
    opt.verbose = 0;
    tstart = tic;
    for trial = 1:maxTrials
        
        %Run EM-BiGAMP  
        [estFinTemp,~,~,estHistEMtemp] = ...
            EMBiGAMP_DL(Q*Y,problem,opt);
        
        %Correct the dictionary
        estFinTemp.Ahat = Q \ estFinTemp.Ahat;
        
        %Reset error function
        opt.error_function = error_function;
        
        %If the sparisty is better and the error is very small or better
        if (sum(sum(estHistEMtemp.p1)) < bestSparsity) && ...
                ( (estHistEMtemp.errZ(end) < bestError)  ||...
                (estHistEMtemp.errZ(end) < -100)  )
            
            %Update
            bestSparsity = sum(sum(estHistEMtemp.p1));
            bestError = estHistEMtemp.errZ(end);
            AhatOptEM = estFinTemp.Ahat;
            estHistEM = estHistEMtemp;
            p1EM = estHistEMtemp.p1;
            
            %Notify user
            %disp(['Accepting new result. Error: ' num2str(bestError)...
            %    '  Average sparsity: ' num2str(bestSparsity/L)...
            %    '  Max Sparsity: ' num2str(max(sum(p1EM)))])
        end
        %keyboard
    end
    tEMGAMP = toc(tstart);
    
    %Save results
    loc = length(results) + 1;
    results{loc}.name = 'EM-BiG-AMP'; %#ok<*AGROW>
    results{loc}.err = estHistEM.errZ(end);
    results{loc}.time = tEMGAMP;
    results{loc}.errHist = estHistEM.errZ;
    results{loc}.timeHist = estHistEM.timing;
    results{loc}.dict = AhatOptEM;
    results{loc}.dictError = dictionary_error_function(results{loc}.dict);
    results{loc}.codingError =...
        coding_error_function(results{loc}.dict);
end




%% SPAMS


if trySpams
    
    %Build params
    spams_param = [];
    spams_param.K = N;
    spams_param.mode = 2;
    spams_param.lambda = 0.1/sqrt(N);
    spams_param.iter = 1000;
    spams_param.verbose = 0;
    
    
    %Trials
    bestSPAMSerror = inf;
    
    tstart = tic;
    for trial = 1:maxTrials
        
        %Do it   
        A_spamsTemp = mexTrainDL(Y,spams_param);
        
        
        SPAMSerrorTemp = dictionary_error_function(A_spamsTemp);
        
        %Update the estimate if this result was superior
        if SPAMSerrorTemp < bestSPAMSerror
            SPAMSerror = SPAMSerrorTemp;
            bestSPAMSerror = SPAMSerror;
            A_spams = A_spamsTemp;
            disp(['Updating Solution. SPAMS Error was: '...
                num2str(SPAMSerror) ' dB'])
        end
        
    end
    tspams = toc(tstart);
    
    
    %Save results
    loc = length(results) + 1;
    results{loc}.name = 'SPAMS'; %#ok<*AGROW>
    results{loc}.err = dictionary_error_function(A_spams);
    results{loc}.time = tspams;
    results{loc}.errHist = results{loc}.err;
    results{loc}.timeHist = zeros(size(results{loc}.errHist));
    results{loc}.dict = A_spams;
    results{loc}.dictError = dictionary_error_function(results{loc}.dict);
    results{loc}.codingError =...
        coding_error_function(results{loc}.dict);
    
    
end


%% DL BCD1


if tryBCD1
    %Trials
    tstart = tic;
    A_bcd = DL_BCD_global(Y, maxTrials, 'L1', .5);
    BCDerror = dictionary_error_function(A_bcd);
    tbcd = toc(tstart);
    disp(['Updating Solution. BCD \tau=.5 Error was: '...
          num2str(BCDerror) ' dB'])
    %Save results
    loc = length(results) + 1;
    results{loc}.name = 'BCD \tau = .5'; %#ok<*AGROW>
    results{loc}.err = dictionary_error_function(A_bcd);
    results{loc}.time = tbcd;
    results{loc}.errHist = results{loc}.err;
    results{loc}.timeHist = zeros(size(results{loc}.errHist));
    results{loc}.dict = A_bcd;
    results{loc}.dictError = dictionary_error_function(results{loc}.dict);
    results{loc}.codingError =...
        coding_error_function(results{loc}.dict);
    
    
end

%% DL BCD2


if tryBCD2
    tstart = tic;
    A_bcd2 = DL_BCD_global(Y, maxTrials, 'rand', inf);
    BCD2error = dictionary_error_function(A_bcd2);
    disp(['Updating Solution. BCD infty Error was: '...
          num2str(BCD2error) ' dB'])
    tbcd2 = toc(tstart);
    
    %Save results
    loc = length(results) + 1;
    results{loc}.name = 'BCD \tau = infinity'; %#ok<*AGROW>
    results{loc}.err = dictionary_error_function(A_bcd2);
    results{loc}.time = tbcd2;
    results{loc}.dict = A_bcd2;
    results{loc}.dictError = dictionary_error_function(results{loc}.dict);
    results{loc}.codingError =...
        coding_error_function(results{loc}.dict);
end


%% ER-SpUD


if tryErspud
    %Call it
    tic
    %[Aerspud,Xerspud] = ER_SpUD_SC(Y);
    [Aerspud,Xerspud]=dl_spud(Y);   %#ok<NASGU>
    tErspud = toc;
    
    %Save results
    loc = length(results) + 1;
    results{loc}.name = 'ER-SpUD (proj)'; %#ok<*AGROW>
    results{loc}.err = dictionary_error_function(Aerspud);
    results{loc}.time = tErspud;
    results{loc}.errHist = results{loc}.err;
    results{loc}.timeHist = zeros(size(results{loc}.errHist));
    results{loc}.dict = Aerspud;
    results{loc}.dictError = dictionary_error_function(results{loc}.dict);
    results{loc}.codingError =...
        coding_error_function(results{loc}.dict);
end


%% K-SVD

if tryKsvd
    
    %Build params
    ksvd_params = [];
    ksvd_params.data = Y;
    ksvd_params.Tdata = K;
    ksvd_params.dictsize = N;
    ksvd_params.iternum = 100;
    ksvd_params.exacty = 1;
    ksvd_params.verbose = 0;
    
    %Trials
    bestKSVDerror = inf;
    tstart = tic;
    for trial = 1:maxTrials
        
        %Do it
        
        [A_ksvdtemp,~,err_ksvdtemp] = ksvd(ksvd_params);
        
        
        %Fix error
        err_ksvdtemp = err_ksvdtemp.^2 * numel(Y);
        
        %Update the estimate if this result was superior
        if err_ksvdtemp(end) < bestKSVDerror
            bestKSVDerror = err_ksvdtemp(end);
            err_ksvd = err_ksvdtemp;
            A_ksvd = A_ksvdtemp;
            disp(['Updating Solution. KSVD Error was: '...
                num2str(10*log10(err_ksvd(end)/norm(Z,'fro')^2))])
        end
        
    end
    tksvd = toc(tstart);
    
    %Save results
    loc = length(results) + 1;
    results{loc}.name = 'K-SVD'; %#ok<*AGROW>
    results{loc}.err = 10*log10(err_ksvd(end)/norm(Z,'fro')^2);
    results{loc}.time = tksvd;
    results{loc}.errHist = 10*log10(err_ksvd/norm(Z,'fro')^2);
    results{loc}.timeHist = zeros(size(err_ksvd));
    results{loc}.dict = A_ksvd;
    results{loc}.dictError = dictionary_error_function(results{loc}.dict);
    results{loc}.codingError =...
        coding_error_function(results{loc}.dict);
    
    
end

%% Store the options structures in results
results{1}.optIn = optIn;
results{1}.trueDict = A;
results{1}.trueEncoding = X;


%% Show Results
if false
disp(['.........................mu:' num2str(optIn.mu) ',K:'...
        num2str(optIn.K) ' SPAMS: ' num2str(SPAMSerror) ' dB'...
        ',BCD:' num2str(BCDerror) ' dB'...
        ' AMP:' num2str(dictionary_error_function(results{1}.dict)) ' dB.'...
        ' K-SVD:' num2str(dictionary_error_function(A_ksvd)) ' dB.'...
        'Erspud:' num2str(dictionary_error_function(Aerspud)) ' dB.'])
end
if nargin == 0
    results{:} %#ok<NOPRT>
    disp('Note that dictError is the normalized error in dB for recovering')
    disp('the dictionary for each algorithm')
end

