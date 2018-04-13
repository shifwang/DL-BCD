len_k = length(K_list);
len_m = length(M_list);
%len_k = 20;
%len_m = 10;
len_trial = 1;
amp_time = nan(len_k, len_m, len_trial);
spams_time = nan(len_k, len_m, len_trial);
bcd1_time = nan(len_k, len_m, len_trial);
bcd2_time = nan(len_k, len_m, len_trial);
spud_time = nan(len_k, len_m, len_trial);
ksvd_error = nan(len_k, len_m, len_trial);
cd sim4_output/
files=dir(fullfile('./','sim4*.mat'));
%addpath('./laplace_sim')
trial_index = 1;
for ind = 1:length(files)
    file = files(ind);
    load(file.name)
    assert(len_k == length(K_list), sprintf('len_k %d not equal to length(K_list) %d ', len_k, length(K_list)))
    assert(len_m == length(M_list), sprintf('len_m %d not equal to length(M_list) %d ', len_m, length(M_list)))
    for i = 1:len_k
        for j = 1:len_m
            try
                if K_list(i) < M_list(j)
                    amp_time(i,j, trial_index) = y{i,j}{1}.time;
                    spams_time(i,j, trial_index) = y{i,j}{2}.time;
                    bcd1_time(i,j, trial_index) = y{i,j}{3}.time;
                    bcd2_time(i,j, trial_index) = y{i,j}{4}.time;
                    spud_time(i,j, trial_index) = y{i,j}{5}.time;
                    ksvd_time(i,j, trial_index) = y{i,j}{6}.time;
                else
                    amp_time(i,j, trial_index)   = 0;
                    spams_time(i,j, trial_index) = 0;
                    bcd1_time(i,j, trial_index)   = 0;
                    bcd2_time(i,j, trial_index)   = 0;
                    spud_time(i,j, trial_index)  = 0;
                    ksvd_error(i,j, trial_index)  = 0;
                end
            catch
            end
        end
    end
    trial_index = trial_index + 1;
end
cd ..
%assert(trial_index == 21)
%%
out = [f(amp_time(10,19,:), 3),
f(spams_time(10,19,:), 3),
f(bcd1_time(10,19,:), 3),
%f(bcd2_time(10,19,:), 3),
f(spud_time(10,19,:), 3),
f(ksvd_time(10,19,:), 3)];
bar(log10(out));
%semilogy()
set(gca,'XTickLabel',{'AMP', 'SPAMS', 'DL-BCD', 'ER-SpUD', 'KSVD'});
ylabel('Running time (sec)')
xlabel('Algorithms')
%% plot
f = @(x,t) nanmean(x,t);
p1 = subplot(2,3,1);
phase_plot(p1, f(amp_time(:,:,:),3), 'EM-BiG-AMP',...
    '$dimension$', 1:length(M_list), M_list,...
    'Number of Nonzeros', 1:length(K_list), K_list)

p2 = subplot(2,3,2);
%p2 = axes;
phase_plot(p2, f(spams_time(:,:,:),3), 'SPAMS',...
    '$dimension$', 1:length(M_list), M_list,...
    'Number of Nonzeros', 1:length(K_list), K_list)
p3 = subplot(2,3,3);
phase_plot(p3, f(bcd1_time(:,:,:),3), 'BCD \tau=.5',...
    '$dimension$',  1:length(M_list), M_list,.../sim4_output/
    'Number of Nonzeros', 1:length(K_list), K_list)
p4 = subplot(2,3,4);
phase_plot(p4, f(spud_time(:,:,:),3), 'ER-SpUD',...
    '$dimension$', 1:length(M_list), M_list,...
    'Number of Nonzeros', 1:length(K_list), K_list)
p5 = subplot(2,3,5);
phase_plot(p5, f(ksvd_error(:,:,:),3), 'K-SVD',...
    '$dimension$', 1:length(M_list), M_list,...
    'Number of Nonzeros', 1:length(K_list), K_list)
p6 = subplot(2,3,6);
phase_plot(p6, f(bcd2_time(:,:,:),3), 'BCD \tau=\infty',...
    '$dimension$',  1:length(M_list), M_list,...
    'Number of Nonzeros', 1:length(K_list), K_list)
