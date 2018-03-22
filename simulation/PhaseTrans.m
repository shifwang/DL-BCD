% Get the phase transition plot
% Reference dictionary is constant linearity with coherence mu
% use dual.m to do L1 minimization
function PhaseTrans(K, samplesize, len, n_trial, storage_path)
    %disp(str2num(getenv('SLURM_CPUS_PER_TASK')))
    %feature('numThreads', 4);
    %addpath(genpath('/accounts/gen/vis/wang.yu/github/DictionaryLearning/matlab/algorithm'))
    addpath(genpath('../algorithm'))
    coherence = linspace(0,.99,len);
    sparsity = linspace(.1,.99,len);
    success_rate = zeros(len, len, n_trial);
    tic;
    for i = 1:len
        for j = 1:len
            for k = 1:n_trial
    	    % set seed
    	    rng(997*i + 887*j + 373*k)
                % generate data
                mu         = coherence(i);
                p          = sparsity(j);
                true_dict  = ((1 - mu) * eye(K) + mu * ones(K))^.5;
                true_coef  = abs(randn(K, samplesize) .* (rand(K, samplesize) < p));
                % run simulation
                [max_dist, ~] = is_sharp(true_coef, true_dict, 1e-2, 1e-5, inf, 0);
                
                success_rate(i, j, k) = log10(max_dist);
            end
        end
    end
    toc;
    % plot
    pp = pcolor(sparsity, coherence, mean(success_rate, 3));
    colormap gray
    set(gca, 'XTickLabel', round(sparsity * 20)/20)
    set(gca, 'YTickLabel', round(coherence*20)/20)
    set(pp, 'EdgeColor', 'none');
    xlabel('p')
    ylabel('\mu')
    hold on;
    x = sparsity;
    y = 2 * x - 1;
    h=plot(x, y, '-');
    legend(h, 'Phase boundary')
    save(storage_path, 'K', 'samplesize', 'len', 'n_trial', 'success_rate', 'sparsity', 'coherence')
end