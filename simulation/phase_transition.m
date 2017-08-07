% script to get theoretical phase transition boundary
addpath('../algorithm/')
len_k = 8;
len_mu = 19;
len_dist = 2;
K = 10;
s_list = 2:(1+len_k);
mu_list = linspace(0.1, 0.9, len_mu);
len_trial = 1;
max_dist = nan(len_k, len_mu, len_trial, len_dist);
sharp_bool = nan(len_k, len_mu, len_dist);
verbose = 1;
samplesize = 2000;
dist_type1 = 1;
thres2 = 1;
test_points = abs(randn(40000,1));
scale = mean((test_points < thres2).*test_points)/sqrt(2/pi);        
for i = 1:len_k
    if verbose
        disp(i);
    end
    for j = 1:len_mu
        s = s_list(i);
        mu = mu_list(j);
        true_dict  = ((1 - mu) * eye(K) + mu * ones(K))^.5;
        for trial_index = 1:len_trial
            true_coef  = randn(K, samplesize);
            for ind = 1:samplesize
                true_coef(randperm(K, K - s),ind) = 0;
            end
            %data       = true_dict * true_coef;
            sharp_bool(i, j, dist_type1) = (sqrt(s) * (mu + 0.01) * scale < (K - s)/(K - 1));
            max_dist(i, j, trial_index, dist_type1) = is_sharp(true_coef, ...
                true_dict, ... 
                1e-2, ... %slope
                1e-6, ... %thres1
                thres2, ... %thres2
                0); %verbose
        end
    end
end
%%
f = @(x,t) nanmean(x,t);
p1 = subplot(2,3,1);
phase_plot(p1, log10(f(max_dist(:,:,:,dist_type1),3)), 'Phase Transition',...
    '$\mu$',2:2:19, .1:.1:.9,...
    'Sparsity', 1:8, 2:9, [-15,3])
p2 = subplot(2,3,4);
phase_plot(p2, -sharp_bool(:,:,dist_type1), 'Phase Transition',...
    '$\mu$',2:2:19, .1:.1:.9,...
    'Sparsity', 1:8, 2:9, [-1, 0])


