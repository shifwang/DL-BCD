d = 10;
N = 100;
p = linspace(0, 1, N);
mu = linspace(0, 1, N);
GS_pred = (1 - p)./sqrt(d - 1);
Wu_pred = zeros(size(p));
Laplace_pred = zeros(size(p));
for i = 1:length(p)
    for s = 0:(d - 1)
        Wu_pred(i) = Wu_pred(i) + nchoosek(d - 1, s)*p(i)^s * (1 - p(i))^(d - s) * sqrt(s);
    end
    Wu_pred(i) = Wu_pred(i)/(p(i) * (d - 1));
    for s = 0:(d - 1)
        Laplace_pred(i) = Laplace_pred(i) + nchoosek(d - 1, s)*p(i)^s * (1 - p(i))^(d - s) * mean(abs(sum(log(rand(s, 1000)./rand(s, 1000)), 1)));
    end
    Laplace_pred(i) = Laplace_pred(i)/(p(i) * (d - 1) * mean(abs(log(rand(1, 1000)./rand(1, 1000)))));
end
h2 = plot(mu, GS_pred, '-r', mu, Wu_pred, '-b', mu, Laplace_pred, '.r', 'LineWidth', 4);
legend(h2, 'Gribonval and Schnass')%, 'Wu and Yu')
xlabel('p')
ylabel('\mu')
title('d = 10, infinite sample size')

%% add empirical plot
addpath(genpath('../'))
local_bool = zeros(N/10, N/10, 5);
for i = 1:N/10
    for j = 1:N/10
        %true_dict  = ((1 - mu(i*10 - 9)) * eye(d) + mu(i*10-9) * ones(d))^.5;
        true_dict = normalize_columns(randn(d, d));
        for trial_index = 1:5
            %true_coef = randn(d, 1000);
            true_coef  = log(rand(d, 1000) ./ rand(d, 1000));
            %true_coef = (randi(2, d, 1000) - 1.5)*2;
            true_coef  = (true_coef .* (rand(d, 1000) < p(j*10-1)));
            local_bool(i, j, trial_index) = is_sharp(true_coef, ...
                true_dict, ... 
                0, ... %slope
                1e-6, ... %thres1
                inf, ... %thres2
                0) < 1e-6; %verbose
        end
    end
end
local_bool = mean(local_bool, 3);
%%
h = pcolor(mu((1:N/10) * 10 - 9), p((1:N/10) * 10 - 9), 1 - local_bool);
colormap gray
set(h, 'EdgeColor', 'none');
axis([0, 1, 0, 1])
hold on;
%xlabel('p');
%ylabel('\mu');
%xticks(mu((1:N/10) * 10 - 5));
%xticklabels(strread(num2str(mu((1:N/10) * 10 - 5)),'%s'));
%yticks(mu((1:N/10) * 10 - 5));
%yticklabels(strread(num2str(mu((1:N/10) * 10 - 5)),'%s'));

