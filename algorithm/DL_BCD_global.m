function [dict, coef, out] = DL_BCD_global(data, MAXTRIAL, init, thres2)
% Complete Dictionary Learning through L1 minimization
% Input:
%   data    : K by N matrix
%   MAXTRIAL: double, maximum number of trials
%   init    : string, how to select initial dict
%   thres   : double, truncated L1 threshold
%
% Output:
%   dict    : K by K matrix, learned dictionary
%   coef    : K by N matrix, coefficient matrix
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
for iter = 1:MAXTRIAL
    if strcmp(init, 'orth')
        options.dict = orth(randn(K,K)); % initial dict
    elseif strcmp(init, 'rand')
        options.dict = normalize_columns(data(:,randperm(N, K)));%FIXME: might not have full rank.
    elseif strcmp(init, 'L1')
        options.dict = orth(randn(K,K)); % initial dict via orthogonal matrix
        options.thres2 = inf;
        [dict, coef, summary] = DL_BCD(data, options);
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
    if summary.max_dist < 1e-9
        break;
    end
end
out.success = summary.max_dist;
out.timing = toc(start);
out.iter = iter;
