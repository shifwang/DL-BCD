function [best_dict, best_coef, out] = DL_BCD_global(data, MAXTRIAL, init, thres2)
% Complete Dictionary Learning through L1 minimization
% Input:
%   data    : K by N matrix
%   MAXTRIAL: double, maximum number of trials
%   init    : string, how to select initial dict
%   thres2  : double, truncated L1 threshold
%
% Output:
%   best_dict    : K by K matrix, learned dictionary
%   best_coef    : K by N matrix, coefficient matrix
%   out          : struct, info of the algorithm
%                    iter: int, how many iterations
%                    success: double, how much distance when perturbing
%                      the dictionary
%                    timing: double, the time of the process.
% by Yu Wang, wang.yu@berkeley.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Configuration  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[K, N] = size(data);
options.MAXITER = 500; % Maximum number of iterations
options.thres1 = 1e-5; % threshold 1 that smooths zero
options.is_sharp = true; % test if the result is sharp
options.method = 'bfgs'; % use bfgs to solve sub-problem
options.verbose = 0; % suppress output
options.random = 0;
options.subsample = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Prepare the output    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out = struct();
out.timing = nan;
out.iter = nan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Main body        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start = tic;
best_dict = 1;
best_coef = 1;
best_max_dist = inf;
for iter = 1:MAXTRIAL
    if strcmp(init, 'orth')
        options.dict = orth(randn(K,K)); % initial dict
    elseif strcmp(init, 'rand')
        options.dict = normalize_columns(randn(K,K));
    elseif strcmp(init, 'randcol')
        options.dict = normalize_columns(data(:,randperm(N, K)));%FIXME: might not have full rank.
    elseif strcmp(init, 'L1')
        options.dict = orth(randn(K,K)); % initial dict via orthogonal matrix
        options.thres2 = inf;
        %options.subsample = .5;
        [dict, ~, ~] = DL_BCD(data, options);
        %options.subsample = 1;
        options.dict = dict;
    elseif strcmp(init, 'spams')
        spams_param = [];
        spams_param.K = N;
        spams_param.mode = 2;
        spams_param.lambda = 0.1/sqrt(N);
        spams_param.iter = 1000;
        spams_param.verbose = 0;
        options.dict = mexTrainDL(data(:,randperm(N)),spams_param);
    else
        error('wrong init value: %s', init);
    end
    options.thres2 = thres2; % threshold 2 that truncates 
    [dict, coef, summary] = DL_BCD(data, options);
    if best_max_dist > summary.max_dist
        best_max_dist = summary.max_dist;
        best_dict = dict;
        best_coef = coef;
    end
    if summary.max_dist < 1e-9
        break;
    end
end
out.success = best_max_dist;
out.timing = toc(start);
out.iter = iter;

