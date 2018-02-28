function [dict, coef, summary] = DL_BCD(data, options)
% Complete Dictionary Learning through L1 minimization
% Input:
%   data    : K by N matrix
%   options : struct, parameter settings.
%             MAXITER - maximum number of iterations
%             dict    - initial dict.
%             thres1  - double, default 1e-5, threshold for recognizing zeros
%             thres2  - double, default inf, threshold to be neglected
%             method  - string, default bfgs, method to optimize the subproblem.
%                           can also be grad, gradient descent.
%             is_sharp- bool, default true, check whether the local minimum is sharp.
%             random  - bool, default false, whether randomly permute
%                 coordinats for each iteration.
%             subsample - double between 0 and 1,  do subsampling
%                 for each iteration, 1 means no subsampling.
%             verbose - int, default 0, the higher the more detailed output
%
% Output:
%   dict    : K by K matrix, learned dictionary
%   coef    : K by N matrix, coefficient matrix
%   summary : struct, information of the algorithm
%             verbose    : bool, whether to output intermediate results
%             obj_val    : 1 by MAXITER matrix, objective for each iteration
%             n_zero     : 1 by MAXITER matrix, percentage of zeros for each iteration
%             flag       : K by MAXITER matrix, exit flag for each iteration
%             timing     : double, time of the whole algorithm
%             perturb    : double, the perturbation made on m to test sharpness
%             max_dist   : double, the L2 distance of the solution when m is perturbed
%             test_timing: double, the timing of the is_sharp test
% by Yu Wang, wang.yu@berkeley.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Configuration  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[K, N]  = size(data);
if isfield(options, 'MAXITER')
    MAXITER = options.MAXITER;
else
    MAXITER = 500;
end
if isfield(options, 'dict')
    dict    = options.dict;
else
    dict = orth(randn(K,K));
end
if isfield(options, 'thres1')
    thres1 = options.thres1;
else
    thres1 = 1e-5;
end
if isfield(options, 'thres2')
    thres2 = options.thres2;
else
    thres2 = inf;
end
if isfield(options, 'is_sharp')
else
    options.is_sharp = true;
end
if isfield(options, 'random')
else
    options.random = false;
end
if isfield(options, 'subsample')
else
    options.subsample = 1;
end
if isfield(options, 'method')
    if strcmp(options.method, 'bfgs')
        bfgs = true;
        gd = false;
    elseif strcmp(options.method, 'grad')
        gd = true;
        bfgs = false;
    else
        error('options.method can only be bfgs or grad, yet received %s', options.method)
    end
else
    bfgs = true;
    gd = false;
end
if isfield(options, 'verbose')
    verbose = options.verbose;
else
    verbose = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Prepare the output    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
summary = struct();
summary.obj_val = nan(1, MAXITER);
summary.n_zero = nan(1, MAXITER);
summary.flag = nan(K, MAXITER);
summary.timing = nan;
summary.verbose = verbose;
summary.perturb = nan;
summary.max_dist = nan;
summary.test_timing = nan;

dual = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Main Body        %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coef    = dict\data;
invD = inv(dict);
start = tic;
for iter = 1:MAXITER
    if verbose > 0
        fprintf('Iter %d', iter)
    end
    if options.random
        ordered = [1, randperm(K-1)+1] ;% always start from 1 to avoid prev ordered[K] = new ordered[1]
    end
    for j = 1:K
        if options.random
            j = ordered(j);
        end
        if options.subsample ~= 1
            ratio = options.subsample + (1 - options.subsample)/MAXITER * iter;
            coef = coef(:,randi(N, 1, floor(N * ratio)));
        end
        if dual
            dual_norms = sqrt(sum(invD.^2));
        end
        % Preparation
        if iter == 1
            %add a small noise around initial dict
            %  to avoid numerical instability.
            x = randn(1, K)*0.01;
        else
            x = zeros(1, K);
        end
        if dual
            x(j) = dual_norms(j);
        else
            x(j) = 1;
        end
        m = dict(:, j)' * dict;
        if max(abs(m)) > 1.01
            warning('m is wrong, is dictionary scaled properly?')
        end
        m(j) = 0;
        if dual
            obj_per_feature = sum(abs(coef).*(abs(coef) < thres2), 2)/N .* dual_norms';
            obj_per_feature(j) = obj_per_feature(j) / dual_norms(j);
        else
            obj_per_feature = sum(abs(coef).*(abs(coef) < thres2), 2)/N;
        end
        %Update dual vector
        % Formulation for j = 1: minimize
        %   E| c1 + c2 x2 + ... + cK xK| + E|c2| sqrt((x2 - m2)^2 + 1 -
        %   m^2) + ... + E|cK| sqrt((xK - mK)^2 + 1 - mK^2)
        % where ci is the coefficients, mi is the cosine of atom i and 1.
        if bfgs
            B = eye(K);
            stepsize = 1;
            ABSTOL   = 1e-14;
            bigger = 1.2;
            smaller = .1;
            percent = 0.1;
            not_one = abs(coef(j,:)) < thres2;
            flag = nan;
            for inner_iter = 1:1000
                if verbose > 1
                    if mod(inner_iter, 100) == 0
                        fprintf('feature %d, inner iter : %d', j, inner_iter);
                    end
                end
                tmp = x * coef;
                iszero = abs(tmp) < thres1;
                aprox_sign = (iszero .* (tmp/thres1) + (~iszero) .* (sign(tmp).*not_one));
                grad1 = aprox_sign * coef' / N;
                grad2 = (x - m) ./ sqrt((x - m).^2 + 1 - m.^2) .* obj_per_feature';
                grad2(j) = -grad1(j); % avoid singularity
                grad = grad1 + grad2;
                if any(~isreal(grad(1)))
                    error('gradient should not be real numbers but complex numbers obtained.');
                end
                if inner_iter > 1
                    yk = grad - old_grad;
                    sk = x - xold;
                    if yk * sk' > 0
                        AA = yk' * yk / (yk * sk') - (B * sk') * sk * B'/ (sk * B * sk');
                    else
                        if verbose > 2
                            warning('yk*sk is negative');
                        end
                        AA = 0;
                    end
                    B = B + AA;
                end
                % Update xold and old_grad 
                xold = x;
                old_grad = grad;
                new_grad = grad/B;
                if isnan(new_grad(1))
                    if verbose > 2
                        warning('new_grad has nan. exiting...');
                    end
                    flag = 6; %has nan, should not happen.
                    break
                end
                bound = new_grad * grad';
                obj_old = mean(abs(tmp.*not_one)) + sqrt((x - m).^2 + 1 - m.^2)*obj_per_feature;
                stepsize = stepsize * bigger;
                xnew = x - stepsize * new_grad;
                obj_new = mean(abs(xnew * coef).*not_one) + sqrt((xnew - m).^2 + 1 - m.^2)*obj_per_feature;
                while obj_new > obj_old - percent * stepsize * bound
                    stepsize = stepsize * smaller;
                    xnew = x - stepsize * new_grad;
                    obj_new = mean(abs(xnew * coef).*not_one) + sqrt((xnew - m).^2 + 1 - m.^2)*obj_per_feature;
                end
                x = xnew;
                if obj_old - obj_new < ABSTOL
                    %disp(obj_old - obj_new)
                    flag = 0; % objective decreases too small
                    break
                end
            end
            if isnan(flag)
                flag = 1; %iteration reaches max, should not happen
            end
        elseif gd
            inner_iter = 0;
            stepsize = 1;
            ABSTOL   = 1e-10;
            bigger = 1.2;
            smaller = .1;
            percent = 0.1;
            not_one = abs(coef(j,:)) < thres2;
            while true                
                if verbose > 1
                    if mod(inner_iter, 100) == 0
                        fprintf('feature %d, inner iter %d', j, inner_iter);
                    end
                end
                tmp = x * coef;
                iszero = abs(tmp) < thres1;
                aprox_sign = (iszero .* (tmp/thres1) + (1 - iszero) .* (sign(tmp).*not_one));
                grad1 = aprox_sign * coef' / N;
                grad2 = (x - m) ./ sqrt((x - m).^2 + 1 - m.^2) .* obj_per_feature';
                grad2(j) = -grad1(j); % avoid singularity
                bound = sum((grad1 + grad2).^2);
                obj_old = mean(abs(tmp.*not_one)) + sqrt((x - m).^2 + 1 - m.^2)*obj_per_feature;
                stepsize = stepsize * bigger;
                xnew = x - stepsize * (grad1 + grad2);
                obj_new = mean(abs(xnew * coef).*not_one) + sqrt((xnew - m).^2 + 1 - m.^2)*obj_per_feature;
                while obj_new > obj_old - percent * stepsize * bound
                    stepsize = stepsize * smaller;
                    xnew = x - stepsize * (grad1 + grad2);
                    obj_new = mean(abs(xnew * coef).*not_one) + sqrt((xnew - m).^2 + 1 - m.^2)*obj_per_feature;
                end
                x = xnew;
                inner_iter = inner_iter + 1;
                if obj_old - obj_new < ABSTOL
                    flag = 0;
                    break
                end
            end
        else
            error('Optimization method not recognized');
        end
        if verbose > 1
            fprintf('stepsize:%6.3e, inner_iter:%d, max gradient:%6.3e',stepsize, inner_iter, max(abs(grad)));
            fprintf('objective %6.3e', obj_new);
        end
        % Update dict and coef
        invD(j,:) = x * invD;
        for h = 1:K
            if h ~= j
                invD(h,:) = invD(h,:) * sqrt((x(h) - m(h))^2 + 1 - m(h)^2);
            end
        end
        if dual
            dict = normalize(invD^-1);
            invD = dict^-1;
        end
        coef = invD * data;
        dict = inv(invD);
        obj  = sum(sum(min(abs(coef), thres2)))/N;
        % summary
        summary.flag(j, iter) = flag;
    end
    if verbose > 0
        fprintf('objective %6.3e', obj);
    end
    summary.obj_val(iter) = obj;
    summary.n_zero(iter) = mean(mean(abs(coef) < thres1));
    if iter > 1
        if ~dual && summary.obj_val(iter-1) - summary.obj_val(iter) < ABSTOL
            break;
        end
    end
end
summary.timing = toc(start);
% check if the obtained dictionary is indeed the reference dictionary
if options.is_sharp
    slope = 1e-2;
    [max_dist, testtime] = is_sharp2(coef, dict, slope, thres1, thres2, verbose);
    summary.perturb = slope;
    summary.max_dist = max_dist;
    summary.test_timing = testtime;
end
if ~isreal(dict(1,1))
    error('dictionary should be real numbers, something is wrong.');
end
%summary
