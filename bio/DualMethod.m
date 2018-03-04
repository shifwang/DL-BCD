function [z, FLAG, iter, betas, obj, err, Lagrangian, timing, err_rel, err_abs] = DualMethod(X, mu, beta, K, num, z0)
  % Dual formulation via ADMM
  %
  % Solves the following problem via ADMM:
  %      minimize    f(X*d) + I(||z||_2 = 1 and z in span(X))
  %      subject to  z = d
  % where f is a relaxed L1.
  %
  % Args:
  %   X         : matrix
  %   mu        : relaxation parameter of L1
  %   beta      : initial penalty parameter in ADMM
  %   K         : K dimension of span(X)
  %   num       : number of vectors learned
  %   z0        : initial value
  %
  % Returns:
  %   z         : estimated vector 
  %   FLAG      : a string indicating the exit status
  %              'Max iteration' : Max iteration has been met.
  %              'Relative error': RELTOL has been met, i.e., 
  %                    ||z^k - d^k||_inf < RELTOL * \| [z^k;d^k] \|_inf
  %              'Absolute error': ABSTOL has been met, i.e.,
  %                    ||z^k - d^k||_inf < ABSTOL
  %              'Unbounded'     : sequence is unbounded, i.e.,
  %                    ||d^k||_inf > LARGE
  %              'beta too large': cannot find appropriate beta, i.e.,
  %                    beta > LARGE (only happens when beta is set automatically)
  %                   
  %   iter       : final iteration
  %   betas      : beta of different iterates
  %   obj        : vector of || X d^k ||_1 for each k
  %   err        : vector of || z^k - d^k ||_inf for each k
  %   Lagrangian : vector of Lagrangian function for each k
  %   timing     : elapsed time of the algorithm
  %   err_rel    : the relative lower bound of different iterates
  %   err_abs    : the absolute lower bound of different iterates
  %
  % More information can be found in the paper linked at:
  %   http://arxiv.org/abs/1511.06324

  % other parameters
  %   x0        : initial point, default is the zero vector
  %   AUTO      : whether beta can be changed by the program
  %   MAXCOUNTS : when Lagrangian increases MAXCOUNTS times, beta = beta * SCALE
  %   SCALE     : beta = beta * SCALE when Lagrangian increases MAXCOUNTS times.
  %   RELTOL    : relative error tolerance 
  %   ABSTOL    : absolute error tolerance 
  %   MAXITER   : max iteration
  %   VERBOSE   : whether print details of the algorithm
  AUTO      = true;
  SCALE     = 1.2;
  RELTOL    = 1e-5;
  ABSTOL    = 1e-6;
  MAXCOUNTS = 5;
  MAXITER   = 1000;
  VERBOSE   = true;
  % Sanity check
  assert(mu > 0);
  assert(beta > 0);
  
  % Default constant
  LARGE  = 1e6;
  %[m, n] = size(A);
  if AUTO
    increasecounts = 0;
  end
  % Main body
  tic; % record the time
  % Initialize
  if ~exist('z0')
      %z0  = X(randi(size(X,1)),:)';
      z0  = randn(size(X,2),1); %Two ways to initialize: arbitrary or use samples
      z0 = z0/norm(z0,2);
  end
  assert(size(X,2) == size(z0,1));
  assert(size(z0,2) == 1);
  z = z0;
  d = z0;
  w = zeros(size(z0));
  obj     = nan(MAXITER, 1);
  err     = nan(MAXITER, 1);
  err_rel = nan(MAXITER, 1);
  err_abs = nan(MAXITER, 1);
  betas   = nan(MAXITER, 1);
  lagrng  = nan(MAXITER, 1);
  % pre-defined functions
  InfNorm    = @(x) max(abs(x));
  f          = @(x) sum(abs(x(abs(x) >= mu))) + sum(x(abs(x)<mu).^2/(2*mu) + mu/2);% relaxed L1
  Lagrangian = @(z,d,w,beta) f(X*d) + w' * (z - d) + beta/2 * sum((z - d).^2);
  if VERBOSE
    fprintf('%4s\t%10s\t%10s\t%10s\n', 'iter', 'obj','lagrng', 'error');
  end
  % Calculate the projection matrix into row space of X
  [U, S, V] = svd(X);
  if ~exist('K')
      S_rank = rank(S);
  else
      S_rank = K;
  end
  P = V(:,1:S_rank) * V(:,1:S_rank)';
  for k = 1:MAXITER 
    % z update
    %------------------------------------------------
    %     Project (d - w/beta) to set {||z||_2 = 1 and z in span(X)}
    %------------------------------------------------
    z = P * (d - w/beta);
    z = z ./ norm(z,2);
    %pause(0.1);
    % d update
    %------------------------------------------------
    %      argmin_d f(X*d) + beta/2 * ||z - d + w/beta||_2^2
    %------------------------------------------------
    d = update_d(d, z + w/beta, X, mu, beta);
    % w update
    w = w + beta * (z - d);
    % record the values
    obj(k)     = f(X*z);
    err(k)     = InfNorm(z - d);
    betas(k)   = beta;
    lagrng(k)  = Lagrangian(z, d, w, beta);
    err_rel(k) = RELTOL * InfNorm([z;d]);
    err_abs(k) = ABSTOL;
    if VERBOSE 
      fprintf('%4d\t%10.4f\t%10.4f\t%10.4f\n', k, obj(k), lagrng(k), err(k));
    end
    % beta update
    if AUTO && k > 1 && lagrng(k) > lagrng(k - 1);
      increasecounts = increasecounts + 1;
    end
    if AUTO && increasecounts > MAXCOUNTS;
      increasecounts = -1;
      beta = beta * SCALE;
    end
    % stopping criteria
    if k == MAXITER
      FLAG = 'Max iteration';
      break;
    end
    if AUTO && beta > LARGE
      FLAG = 'beta too large';
      break;
    end
    if InfNorm(d) > LARGE
      FLAG = 'Unbounded';
      break;
    end
    if err(k) < ABSTOL
      FLAG = 'Absolute error';
      break;
    end
    if err(k) < RELTOL * InfNorm([d;z])
      FLAG = 'Relative error';
      break;
    end
  end
  iter = k;
  timing = toc;
  if VERBOSE
    fprintf('ADMM has stopped at iter %4d because of %10s.\n',k,FLAG);
    fprintf('Elapsed time is %8.1e seconds .\n',timing);
  end
end

function d = update_d(d, z, X, mu, beta)
    % solve the x update
    %    argmin_d f(X*d) + beta/2 * ||d - z||_2^2
    % via Newton's method; for a single subsystem only.
    alpha = 0.1;
    BETA  = 0.5;
    TOLERANCE = 1e-7;
    MAX_ITER = 50;
    f          = @(x) sum(abs(x(abs(x) >= mu))) + sum(x(abs(x)< mu).^2/(2*mu) + mu/2);% relaxed L1
    df         = @(x) sign(x).*(abs(x) >= mu) + (x./mu) .* (abs(x) < mu);
    ddf        = @(x) diag(1/mu.*(abs(x) < mu));
    for iter = 1:MAX_ITER
        fx = f(X*d) + beta/2*sum((d - z).^2);
        g = X'*df(X*d) + beta*(d - z);
        H = X'*ddf(X*d)*X + beta*eye(length(d));
        dx = -H\g;   % Newton step
        dfx = g'*dx; % Newton decrement
        if abs(dfx) < TOLERANCE
            break;
        end
        % backtracking
        t = 1;
        while f(X*(d + t*dx)) + beta/2*sum((d - z + t*dx).^2) > fx + alpha*t*dfx
            t = BETA*t;
        end
        d = d + t*dx;
    end
end
%% Visualization 
%    figure;
%    for i = 1:min(10,12)
%        subplot(4, 3, i)
%        min_value = min(x(:,i));
%        max_value = max(x(:,i));
%        imshow(reshape((x(:,i)-min_value)./(max_value - min_value),32, 32));
%    end

