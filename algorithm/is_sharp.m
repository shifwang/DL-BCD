function [max_dist, time] = is_sharp(coef, dict, slope, thres1, thres2, verbose)
if ~exist('slope')
    slope = 1e-2;
end
test = tic;
max_dist = 0;
K = size(coef, 1);
N = size(coef, 2);
for j = 1:K
    % Preparation
    x = randn(1, K)*1e-3;
    x(j) = 1;
    obj_per_feature = sum(abs(coef).*(abs(coef) < thres2), 2)/N;
    m = dict(:, j)' * dict;
    m = m + slope;
    m(j) = 0;
    %Update dual vector
    % Formulation for j = 1: minimize
    %   E| c1 + c2 x2 + ... + cK xK| + E|c2| sqrt((x2 - m2)^2 + 1 -
    %   m^2) + ... + E|cK| sqrt((xK - mK)^2 + 1 - mK^2)
    % where ci is the coefficients, mi is the cosine of atom i and 1.
    B = eye(K);
    stepsize = 1;
    ABSTOL   = 1e-10;
    bigger = 1.2;
    smaller = .1;
    percent = 0.1;
    not_one = abs(coef(j,:)) < thres2;
    for inner_iter =1:1000
        tmp = x * coef;
        iszero = abs(tmp) < thres1;
        aprox_sign = (iszero .* (tmp/thres1) + (~iszero) .* (sign(tmp).*not_one));
        grad1 = aprox_sign * coef' / N;
        grad2 = (x - m) ./ sqrt((x - m).^2 + 1 - m.^2) .* obj_per_feature';
        grad2(j) = -grad1(j); % avoid singularity
        grad = grad1 + grad2;
        if inner_iter > 1
            yk = grad - old_grad;
            sk = x - xold;
            if yk * sk' > 0
                AA = yk' * yk / (yk * sk') - (B * sk') * sk * B'/ (sk * B * sk');
            else
                %warning('yk*sk is negative')
                AA = 0;
            end
            B = B + AA;
        end
        % Update xold and old_grad
        xold = x;
        old_grad = grad;
        new_grad = grad/B;
        if isnan(new_grad(1))
            %warning('nan');
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
        inner_iter = inner_iter + 1;
        x = xnew;
        if obj_old - obj_new < ABSTOL
            break
        end
    end
    xnew(j) = 0;
    if sum(xnew.^2) > max_dist
        max_dist = sum(xnew.^2);
    end
end
if verbose > 0
    if max_dist <  1e-6
        disp('The obtained dictionary is the sharp local min.')
    else
        disp('The obtained dictionary is NOT the sharp local min.')
    end
end
time = toc(test);
end
