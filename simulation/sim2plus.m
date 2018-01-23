addpath(genpath('../AMP/'))
addpath('../algorithm/')
addpath(genpath('../externs/'))
%%
middles = [];
for K = 10:5:20
N  = 2 * K; % number of candidate sample sizes
N_trials = 40; % number of trials
success1 = nan(N, N_trials);
s = round(.5 * K);
mu = ((K - s)/(K - 1) - 0.2)/sqrt(s);
for i = 1:N    
    samplesize = 20 * K + (i - K)*10;
    for j = 1:N_trials
    true_coef  = randn(K, samplesize);
    for ind = 1:samplesize
        true_coef(randperm(K, K - s),ind) = 0;
    end
    true_dict  = ((1 - mu) * eye(K) + mu * ones(K))^.5;
    %true_dict = normalize(randn(K,K));
    data       = true_dict * true_coef;
    [max_dist, testtime_once] = is_sharp(true_coef,...
        true_dict, ...
        0.01, ... % slope
        1e-6, ... % thres1
        inf, ... %thres2
        0); %verbose
    success1(i, j) = max_dist < 1e-6;
    end
end
f = @(p, x) exp(p(1) * x + p(2))./(1 + exp(p(1)*x + p(2)));
x = linspace(50, N*50, N);
y = mean(success1,2);
para = nlinfit(x', y, f, [.01, -6]);
middles = [middles, round(-para(2)/para(1))];
disp(middles)
end

