d = 20;
N = 100;
p = linspace(0, 1, N);
mu = linspace(0, 1, N);
GS_pred = (1 - p)./sqrt(d - 1);
Wu_pred = zeros(size(p));
Laplace_pred = zeros(size(p));
non_negative_pred = 2 * p - 1;
for i = 1:length(p)
    for s = 0:(d - 1)
        Wu_pred(i) = Wu_pred(i) + nchoosek(d - 1, s)*p(i)^s * (1 - p(i))^(d -1 - s) * sqrt(s);
    end
    Wu_pred(i) = Wu_pred(i)/(p(i) * (d - 1))*(1 - p(i));
    for s = 0:(d - 1)
        %Laplace_pred(i) = Laplace_pred(i) + nchoosek(d - 1, s)*p(i)^s * (1 - p(i))^(d -1 - s) * mean(abs(sum(log(rand(s, 20000)./rand(s, 20000)), 1)));
        %Laplace_pred(i) = Laplace_pred(i) + nchoosek(d - 1, s)*p(i)^s * (1 - p(i))^(d -1 - s) * mean(abs(gamrnd(s, 1, [1, 20000]) - gamrnd(s, 1, [1, 20000])));
        f = @(x, y) abs(x - y).*(x.*y).^(s - 1) .* exp(-x - y)./gamma(s)^2;
        Laplace_pred(i) = Laplace_pred(i) + nchoosek(d - 1, s)*p(i)^s * (1 - p(i))^(d -1 - s) * integral2(f, 0, inf, 0, inf);
    end
    %Laplace_pred(i) = Laplace_pred(i)/(p(i) * (d - 1) * mean(abs(log(rand(1, 20000)./rand(1, 20000)))))* (1 - p(i));
    Laplace_pred(i) = Laplace_pred(i)/(p(i) * (d - 1))* (1 - p(i));
end
%%
subplot(1,2,1)
h2 = plot(mu, Laplace_pred, '-r', mu, Wu_pred, '-b', 'LineWidth', 2);
legend(h2, 'sparse Laplace', 'sparse Gaussian')
%%
subplot(1,2,2)
h3 = plot(mu, Laplace_pred, '-r', mu, Wu_pred, '-b', 'LineWidth', 2);
legend(h3, 'sparse Laplace', 'sparse Gaussian')

%h2 = plot(mu, Laplace_pred, '.-y', 'LineWidth', 4);
%legend(h2, 'new bound')
%h2 = plot(mu, GS_pred, '-r', mu, Wu_pred, '-b', 'LineWidth', 4);%, mu, Laplace_pred, '.-y'
%legend(h2, 'Gribonval and Schnass', 'Wu and Yu')
%h2 = plot(mu, GS_pred, '-r', 'LineWidth', 4)%, mu, Wu_pred, '-b', 'LineWidth', 4);%, mu, Laplace_pred, '.-y'
%legend(h2, 'Gribonval and Schnass')
xlabel('p')
ylabel('\mu')
title('K = 10, n = 1000')

%% add empirical plot
addpath(genpath('../'))
rng(10)
M = 10;
sparse_mu = linspace(0, .95, M);
sparse_p = linspace(0.05, 1, M);
local_bool = zeros(M, M, 5);
for i = 1:M
    disp(i)
    for j = 1:M
        true_dict  = ((1 - sparse_mu(i)) * eye(d) + sparse_mu(i) * ones(d))^.5;
        %true_dict = normalize_columns(randn(d, d));
        for trial_index = 1:5
            %true_coef = randn(d, 1000);
            true_coef = abs(randn(d, 2000));
            %true_coef  = log(rand(d, 1000) ./ rand(d, 1000));
            %true_coef = (randi(2, d, 1000) - 1.5)*2;
            true_coef  = (true_coef .* (rand(d, 2000) < sparse_p(j)));
            options = struct();
            options.thres1 = 1e-6;
            options.thres2 = inf;
            options.is_sharp = 0;
            options.dict = true_dict;
            options.MAXITER = 3;
            [dict, ~, ~] = DL_BCD(true_dict * true_coef, options);
            local_bool(i, j, trial_index) = norm(dict - options.dict, 'F');
%             local_bool(i, j, trial_index) = is_sharp(true_coef, ...
%                 true_dict, ... 
%                 0, ... %slope
%                 1e-6, ... %thres1
%                 inf, ... %thres2
%                 0) < 1e-6; %verbose
        end
    end
end
local_bool = mean(local_bool, 3);
%%
h = pcolor(sparse_mu, sparse_p, log10(local_bool));
colorbar
colormap gray
set(h, 'EdgeColor', 'none');
axis([0, .95, 0.05, 1])
hold on;
%xlabel('p');
%ylabel('\mu');
%xticks(mu((1:N/10) * 10 - 5));
%xticklabels(strread(num2str(mu((1:N/10) * 10 - 5)),'%s'));
%yticks(mu((1:N/10) * 10 - 5));
%yticklabels(strread(num2str(mu((1:N/10) * 10 - 5)),'%s'));

%% sparse Gaussian vs sparse Laplace

