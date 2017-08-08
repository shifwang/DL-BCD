% test how many samples will make sure the reference dictionary is a sharp
% local min
addpath(genpath('../AMP/'))
addpath('../algorithm/')
K  = 20; % dimension
N  = 20; % number of candidate sample sizes
N_trials = 40; % number of trials
success1 = nan(N, N_trials);
%samplesize = 2000;        
s = 10;
mu = .15;
for i = 1:N
    samplesize = i*50+1500;
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
%%
subplot(1,2,1)
%plot(linspace(50, N*50, N), mean(success1,2))
x = linspace(50, N*50, N);
y = mean(success1,2);
%p = [[x']\y;0]';
[fity, ci] = binofit(sum(success1,2), N_trials, .99);
%errorbar(x, fity, ...
%    ci(:,1) - fity, ci(:,2) - fity,...
%    'LineStyle','None',...
%    'Marker', '.', 'MarkerSize', 15)
axis([0, 1010, 0, 1])
hold on;
plot(x, fity, '-', 'LineWidth', 2);
%legend({'timing (\mu \pm \sigma)', 'linear reg'}, ...
%    'location', 'NorthWest', ...
%    'FontSize', 15)
ylabel('Percentage')
xlabel('Samplesize')

