%len_k = length(K_list);
%len_mu = length(mu_list);
len_k = 10;
len_m = 19;
len_trial = 1;
amp_error = nan(len_k, len_m, len_trial);
spams_error = nan(len_k, len_m, len_trial);
bcd_error = nan(len_k, len_m, len_trial);
spud_error = nan(len_k, len_m, len_trial);
ksvd_error = nan(len_k, len_m, len_trial);
files=dir(fullfile('./','sim4*.mat'));
%addpath('./laplace_sim')
trial_index = 1;
for ind = 1:length(files)
    file = files(ind);
    load(file.name)
    assert(len_k == length(K_list))
    assert(len_m == length(mu_list))
    for i = 1:len_k
        for j = 1:len_m
            try
                amp_error(i,j, trial_index) = y{i,j}{1}.dictError;
                spams_error(i,j, trial_index) = y{i,j}{1}.dictError;
                bcd_error(i,j, trial_index) = y{i,j}{2}.dictError;
                spud_error(i,j, trial_index) = y{i,j}{3}.dictError;
                ksvd_error(i,j, trial_index) = y{i,j}{3}.dictError;
            catch
            end
        end
    end
    trial_index = trial_index + 1;
end
assert(trial_index == 21)
%% plot
f = @(x,t) nanmean(x,t);
p1 = subplot(2,3,1);
phase_plot(p1, f(amp_error(:,:,:),3), 'EM-BiG-AMP',...
    '$\mu$',2:2:19, .1:.1:.9,...
    'Sparsity', 1:8, 2:9)

p2 = subplot(2,3,2);
%p2 = axes;
phase_plot(p2, f(spams_error(:,:,:),3), 'SPAMS',...
    '$\mu$',2:2:19, .1:.1:.9,...
    'Sparsity', 1:8, 2:9)
p3 = subplot(2,3,3);
phase_plot(p3, f(bcd_error(:,:,:),3), 'BCD',...
    '$\mu$',2:2:19, .1:.1:.9,...
    'Sparsity', 1:8, 2:9)
p4 = subplot(2,3,4);
phase_plot(p4, f(spud_error(:,:,:),3), 'ER-SpUD',...
    '$\mu$',2:2:19, .1:.1:.9,...
    'Sparsity', 1:8, 2:9)
p5 = subplot(2,3,5);
phase_plot(p5, f(ksvd_error(:,:,:),3), 'K-SVD',...
    '$\mu$',2:2:19, .1:.1:.9,...
    'Sparsity', 1:8, 2:9)