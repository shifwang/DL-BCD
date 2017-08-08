% Simulation 1:
%   Test empirical time complexity of Algorithm 1
addpath(genpath('../AMP/'))
addpath('../algorithm/')
%% how test time scale with N
K  = 20; % dimension
N  = 10; % number of candidate sample sizes
N_trials = 20; % number of trials
dist = nan(N, N_trials);
testtime1 = nan(N, N_trials);
%samplesize = 2000;        
p = 0.7;
mu = 0.5;
for i = 1:N
    samplesize = i*500;
    for j = 1:N_trials
    true_coef  = randn(K, samplesize) .* (rand(K, samplesize) < p);
    true_dict  = ((1 - mu) * eye(K) + mu * ones(K))^.5;
    %true_dict = normalize(randn(K,K));
    data       = true_dict * true_coef;
    [max_dist, testtime_once] = is_sharp(true_coef,...
        true_dict, ...
        0.01, ... % slope
        1e-6, ... % thres1
        inf, ... %thres2
        0); %verbose
    dist(i, j) = max_dist;
    testtime1(i, j) = testtime_once;
    end
end
%%
subplot(1,2,1)
plot(linspace(500, N*500, N), mean(testtime1,2))
x = linspace(500, N*500, N);
y = mean(testtime1,2);
%p = [[x']\y;0]';
p = polyfit(x, y', 1);
fity = polyval(p, 0:5100);
errorbar(x, mean(testtime1',1), ...
    std(testtime1', 1),...
    'LineStyle','None',...
    'Marker', '.', 'MarkerSize', 15)
axis([0, 5100, 0.15, 0.45])
hold on;
plot(0:5100, fity, 'r--', 'LineWidth', 2);
legend({'timing (\mu \pm \sigma)', 'linear reg'}, ...
    'location', 'NorthWest', ...
    'FontSize', 15)
ylabel('Timing (sec)')
xlabel('Samplesize')

%% how test time scales with K
N  = 10;
N_trials = 20; % number of trials
dist = nan(N, N_trials);
testtime = nan(N, N_trials);
samplesize = 400;
p = 0.7;
for i = 1:N
    for j = 1:N_trials
    K = i * 5;
    true_coef  = randn(K, samplesize) .* (rand(K, samplesize) < p);
    true_dict  = ((1 - mu) * eye(K) + mu * ones(K))^.5;
    %true_dict = normalize(randn(K,K));
    data       = true_dict * true_coef;
        [max_dist, testtime_once] = is_sharp(true_coef,...
        true_dict, ...
        0.01, ... % slope
        1e-6, ... % thres1
        inf, ... %thres2
        0); %verbose
    dist(i, j) = max_dist;
    testtime(i, j) = testtime_once;
    end
end
%%
subplot(1,2,2)
x = linspace(5, N*5, N);
y = mean(testtime,2);
p = [[x'.^2, x']\y;0]';
fity = polyval(p, 0:55);
errorbar(x, mean(testtime',1), ...
    std(testtime', 1),...
    'LineStyle','None',...
    'Marker', '.', 'MarkerSize', 15)
axis([0, 55, 0, 2.5])
hold on;
plot(0:55, fity, 'r--', 'LineWidth', 2);
legend({'Timing (\mu \pm \sigma)', 'quadratic reg'}, ...
    'location', 'NorthWest', ...
    'FontSize', 15)
xlabel('Dimension')
ylabel('Timing (sec)')
