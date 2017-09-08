% test DL_BCD_global performance
addpath(genpath('../AMP/'))
addpath('../algorithm/')
addpath(genpath('../externs/'))
%%
K  = 10; % dimension
N  = 1; % number of candidate sample sizes
N_trials = 10; % number of trials
dist = nan(N, N_trials);
success = nan(N, N_trials);
testtime1 = nan(N, N_trials);
iter1 = nan(N, N_trials);
%samplesize = 2000;        
s = 7;
mu = 0.7;
for i = 1:N
    samplesize = 1000;
    for j = 1:N_trials
        disp(j);
    true_coef  = randn(K, samplesize);
    for ind = 1:samplesize
        true_coef(randperm(K, K - s),ind) = 0;
    end
    true_dict  = ((1 - mu) * eye(K) + mu * ones(K))^.5;
    
    %Dictionary error function
    dictionary_error_function =...
        @(q) 20*log10(norm(true_dict -...
            q * find_permutation(true_dict,q),'fro')/norm(true_dict,'fro'));
    data       = true_dict * true_coef;
    [dict, coef, out] = DL_BCD_global(data, 5, 'orth');
    dist(i, j) = dictionary_error_function(dict)
    testtime1(i, j) = out.timing;
    iter1(i, j) = out.iter;
    success(i, j) = out.success;
    end
end
disp(dist);
disp(iter1);
disp(testtime1);
