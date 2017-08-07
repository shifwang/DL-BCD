%len_k = length(K_list);
%len_mu = length(mu_list);
len_k = 8;
len_mu = 19;
len_trial = 1;
is_pos = 1;
amp_error = nan(len_k, len_mu, len_trial, is_pos);
spams_error = nan(len_k, len_mu, len_trial, is_pos);
bcd_error = nan(len_k, len_mu, len_trial, is_pos);
spud_error = nan(len_k, len_mu, len_trial, is_pos);
ksvd_error = nan(len_k, len_mu, len_trial, is_pos);
files=dir(fullfile('./','*.mat'));
%addpath('./laplace_sim')
trial_index_pos = 1;
trial_index = 1;
for ind = 1:length(files)
    file = files(ind);
    load(file.name)
    assert(len_k == length(K_list))
    assert(len_mu == length(mu_list))
    if optIn.pos == 0
        for i = 1:len_k
            for j = 1:len_mu
                try
                amp_error(i,j, trial_index_pos, 1) = y{i,j}{1}.dictError;
                spams_error(i,j, trial_index_pos, 1) = y{i,j}{1}.dictError;
                bcd_error(i,j, trial_index_pos, 1) = y{i,j}{2}.dictError;
                spud_error(i,j, trial_index_pos, 1) = y{i,j}{2}.dictError;
                ksvd_error(i,j, trial_index_pos, 1) = y{i,j}{2}.dictError;
                catch
                end
            end
        end
        trial_index_pos = trial_index_pos + 1;
    else
        for i = 1:len_k
            for j = 1:len_mu
                amp_error(i,j, trial_index, 2) = y{i,j}{1}.dictError;
                spams_error(i,j, trial_index, 2) = y{i,j}{1}.dictError;
                bcd_error(i,j, trial_index, 2) = y{i,j}{1}.dictError;
                spud_error(i,j, trial_index, 2) = y{i,j}{1}.dictError;
                ksvd_error(i,j, trial_index, 2) = y{i,j}{1}.dictError;
            end
        end
        trial_index = trial_index + 1;
    end
end
assert(trial_index == 21)
assert(trial_index_pos == 21)
%% plot
pos_ind = 1;
f = @(x,t) nanmean(x,t);
p1 = subplot(2,3,1);
phase_plot(p1, f(amp_error(:,:,:,pos_ind),3), 'EM-BiG-AMP',...
    '$\mu$',2:2:19, .1:.1:.9,...
    'Sparsity', 1:8, 2:9)

p2 = subplot(2,3,2);
%p2 = axes;
phase_plot(p2, f(spams_error(:,:,:,pos_ind),3), 'SPAMS',...
    '$\mu$',2:2:19, .1:.1:.9,...
    'Sparsity', 1:8, 2:9)
p3 = subplot(2,3,3);
phase_plot(p3, f(bcd_error(:,:,:,pos_ind),3), 'BCD',...
    '$\mu$',2:2:19, .1:.1:.9,...
    'Sparsity', 1:8, 2:9)
p4 = subplot(2,3,4);
phase_plot(p4, f(spud_error(:,:,:,pos_ind),3), 'ER-SpUD',...
    '$\mu$',2:2:19, .1:.1:.9,...
    'Sparsity', 1:8, 2:9)
p5 = subplot(2,3,5);
phase_plot(p5, f(ksvd_error(:,:,:,pos_ind),3), 'K-SVD',...
    '$\mu$',2:2:19, .1:.1:.9,...
    'Sparsity', 1:8, 2:9)