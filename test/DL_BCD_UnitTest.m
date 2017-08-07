addpath('../algorithm/')
addpath(genpath('../AMP'))
%Dictionary error function
dictionary_error_function =...
    @(A, q) 20*log10(norm(A -...
    q * find_permutation(A,q),'fro')/norm(A,'fro'));

samplesize = 1000;
K          = 2;
p          = 0.9;
true_coef  = randn(K, samplesize) .* (rand(K, samplesize) < p);
theta1     = pi/3;
theta2     = 0;
true_dict  = [cos(theta1), cos(theta2); sin(theta1), sin(theta2)];
data       = true_dict * true_coef;
MAXITER    = 100;
options.MAXITER = MAXITER;
options.thres2 = 10;
options.dict = orth(randn(K,K));
[dict, coef, summary] = DL_BCD(data, options);
plot(1:MAXITER, summary.obj_val, 1:MAXITER, repmat(sum(sum(abs(true_coef)))/samplesize, MAXITER), 'r')
figure()
plot(1:MAXITER, summary.n_zero, 'r')
disp(dict)
disp(true_dict)
disp(summary.timing)

%% Test 1 pair
samplesize = 1000;
K          = 2;
p          = 0.2;
true_coef  = abs(randn(K, samplesize) .* (rand(K, samplesize) < p));
theta1     = pi/3;
theta2     = 0;
true_dict  = [cos(theta1), cos(theta2); sin(theta1), sin(theta2)];
data       = true_dict * true_coef;
MAXITER    = 100;
options.MAXITER = MAXITER;
options.dict = orth(randn(K,K));
options.penalty = 0;
options.thres2 = 100;
[dict, coef, summary] = DL_BCD(data, options);
bias = 1;
plot(1:MAXITER, summary.obj_val, 1:MAXITER, repmat(sum(sum(abs(true_coef)))/samplesize, MAXITER), 'r')
disp(dict)
disp(true_dict)
disp(summary.timing)
%% Test 2
samplesize = 1000;
K          = 10;
p          = 0.5;
true_coef  = randn(K, samplesize) .* (rand(K, samplesize) < p);
mu         = 0.24;
true_dict  = ((1 - mu) * eye(K) + mu * ones(K))^.5;
data       = true_dict * true_coef;
MAXITER    = 60;
options.MAXITER = MAXITER;
options.dict = orth(randn(K,K));
options.thres1 = 1e-6;
options.thres2 = inf;
[dict, coef, summary] = DL_BCD(data, options);
plot(1:MAXITER, summary.obj_val, 1:MAXITER, repmat(sum(sum(abs(true_coef)))/samplesize, MAXITER), 'r')
disp(dictionary_error_function(true_dict, dict))
%disp(dict)
%disp(true_dict)
disp(summary.timing);
%% Test 2 pair
samplesize = 2000;
K          = 10;
p          = 0.3;
true_coef  = randn(K, samplesize) .* (rand(K, samplesize) < p) .* repmat(randn(K, 1), 1, samplesize);
mu         = 0.2;
%true_dict  = ((1 - mu) * eye(K) + mu * ones(K))^.5;
true_dict = normalize(randn(K,K));
data       = true_dict * true_coef;
MAXITER    = 60;
options.MAXITER = MAXITER;
options.dict = orth(randn(K,K));
options.thres2 = inf;
options.thres1 = 1e-6;
tic;
[dict, coef, summary] = DL_BCD(data, options);
toc;
plot(1:MAXITER, summary.obj_val, 1:MAXITER, repmat(sum(sum(abs(true_coef)))/samplesize, MAXITER), 'r')
disp(dictionary_error_function(true_dict, dict))
%disp(dict)
%disp(true_dict)

%% Test case: toy image example
% first run generateTestExample.m
[proj_data, tmp] = pca(data, 'center', false);
proj_data = proj_data(:, 1:dict_size)';
options.MAXITER = 100;
options.dict = orth(randn(dict_size,dict_size));
options.thres2 = inf;
[newdict, newcoef, summary] = DL_BCD(proj_data, options);
disp(summary.timing)
disp(summary.max_dist)
disp(summary.test_timing)
%plot newdict
finaldict = tmp(:,1:dict_size) * newdict;
figure;
for i = 1:dict_size
    min_value = min(finaldict(:,i));
    max_value = max(finaldict(:,i));
    subplot(3, 3, i)
    imshow(reshape((finaldict(:,i)-min_value)./(max_value - min_value),width, height));
end
% compare with lasso type
%% configure the parameters
clear param
param.K         = dict_size;  % dictionary size
param.iter      = 3000;       % number of iterations
param.mode      = 2;          % 1 if basis pursuit; 2 if lasso 
param.lambda    = .05;       % tolerance
param.posAlpha  = 0;          % true if the coefficients are positive
param.posD      = 0;          % true if D is non-negative
param.modeD     = 0;          % 0 if all the atoms' length is no bigger than 1.
param.whiten    = 0;          % I don't know what is that
param.pos       = 0;          % to ensure positive constraints when using mexLasso
param.modeParam = 0;          % optimization parameter, I don't know its meaning
param.verbose   = 0;
%% 
UU=mexTrainDL(proj_data,param);
VV=mexLasso(proj_data,UU,param);
UU = tmp(:,1:dict_size) * UU;
figure;
for i = 1:dict_size
    min_value = min(UU(:,i));
    max_value = max(UU(:,i));
    subplot(3, 3, i)
    imshow(reshape((UU(:,i)-min_value)./(max_value - min_value),width, height));
end



