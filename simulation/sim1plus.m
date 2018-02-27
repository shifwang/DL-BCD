% Simulation 1:
%   Test empirical time complexity of Algorithm 1
addpath(genpath('../AMP/'))
addpath('../algorithm/')
%% sensitivity analysis of perturbation level rho
N  = 1; % number of candidate feature sizes
N_trials = 20; % number of trials
N_rho = 10;
sharp = nan(N_trials, N_rho, N);
non_sharp = nan(N_trials, N_rho, N);
for i = 1:N
  K = i * 20;
  s = round(.5 * K);
  % success
  samplesize = 80 * K;
  mu = ((K - s)/(K - 1) - 0.2)/sqrt(s);
  for j = 1:N_trials
    true_coef  = randn(K, samplesize);
    for ind = 1:samplesize
      true_coef(randperm(K, K - s),ind) = 0;
    end
    true_dict  = ((1 - mu) * eye(K) + mu * ones(K))^.5;
    %true_dict = normalize(randn(K,K));
    data       = true_dict * true_coef;
    for k = 1:N_rho
      [max_dist, testtime_once] = is_sharp2(true_coef,...
        true_dict, ...
        k * .05, ... % slope
        1e-6, ... % thres1
        inf, ... %thres2
        0); %verbose
      sharp(j, k, i) = max_dist;
    end
  end
  mu = ((K - s)/(K - 1) + 0.1)/sqrt(s);
  for j = 1:N_trials
    true_coef  = randn(K, samplesize);
    for ind = 1:samplesize
      true_coef(randperm(K, K - s),ind) = 0;
    end
    true_dict  = ((1 - mu) * eye(K) + mu * ones(K))^.5;
    %true_dict = normalize(randn(K,K));
    data       = true_dict * true_coef;
    for k = 1:N_rho
      [max_dist, testtime_once] = is_sharp2(true_coef,...
        true_dict, ...
        k * .05, ... % slope
        1e-6, ... % thres1
        inf, ... %thres2
        0); %verbose
      non_sharp(j, k, i) = max_dist;
    end
  end
end
%%
%save('sim2plus.mat')
subplot(1,2,1)
boxplot(non_sharp(:,:,1), 'Labels', {.05, .1, .15, .2, .25, .3, .35, .4, .45, .5});
xlabel('\rho')

set(gca,'yscale','log')
ylabel('distance')
ylim([min(min(sharp(:,:,1))), max(max(non_sharp(:,:,1)))])
%hold on
subplot(1,2,2)
boxplot(sharp(:,:,1), 'Labels', {.05, .1, .15, .2, .25, .3, .35, .4, .45, .5});
ylim([min(min(sharp(:,:,1))), max(max(non_sharp(:,:,1)))])
xlabel('\rho')
set(gca,'yscale','log')