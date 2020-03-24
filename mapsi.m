function [ w, wse, in_verts, simp, lambdas, scores, ...
    minlambda, minlambdase, A, wcovse ] = mapsi( data, opts, A)
% MAPSI performs Bayesian estimation of the orientations of particles in
% solution from scattering data
%
% INPUT
% data - a mapsi_data structure
% opts - a mapsi_options structure
% A - the kernel matrix, if it is already calculated
%
% OUTPUT
% w - MAP estimate for the weights
% wse - MAP estimate for the weights with the +1SE rule
% in_verts - matrix of spherical coordinates for simplex vertices
% simp - a matrix where each row contains the indices of in_verts
%   associated with each simplex
% lambdas - vector of lambda values used for cross validation
% scores - a matrix of scores for each lambda and each MCCV fold
% minlambda - the optimal value of lambda from MCCV
% minlambdase - the lambda value corresponding to the +1SE rule
% A - the calculated kernel matrix
% wcovse - the sampled covariance matrix of the +1SE weights

% disable warnings
orig_warnings = warning;
warning('off','all')

if isempty(gcp('nocreate'))
    parpool;
    fprintf('\n')
end

rng('shuffle');

% data queue for logging progress
data_queue = parallel.pool.DataQueue;
afterEach(data_queue, @mapsi_waitbar);
data_queue2 = parallel.pool.DataQueue;
afterEach(data_queue2, @(s) mapsi_waitbar(-1,s));

mapsi_waitbar(0);
mapsi_waitbar(1, 'Preprocessing');

% define constants
M = length(data.I);

% get the kurihara mesh
[in_verts, simp] = create_kurihara_mesh(opts.N);
N = length( in_verts );

% calculate A matrix and scale so basis functions are normalized
if nargin < 3
    A = integrate_basis_mesh_2D(data.kernel, simp, in_verts, data.Q, opts.n_quad, 1);
end

Aeq = integrate_basis_mesh_2D(@(in,out) 1, simp, in_verts, 0, opts.n_quad);
beq = 1;

w0 = 1 * ones(N,1) / sum(Aeq);

% calculate the regularization matrix so that w' * B * w is the integral of the
% square magnitude of the gradient of the distribution on the surface 
B = integrate_square_gradient( in_verts, simp, opts.n_quad );

% ------------------ do cross validation to find regularization ---------
lambdas = opts.lambdas/N;     % range of lambdas to use
k = opts.folds;             % number of folds

mapsi_waitbar(1, 'Processing');
mapsi_waitbar(2, k*length(lambdas) + 2);

inds = crossvalinds(M, k);

scores = zeros(k, length(lambdas));

traininds = cell(k,1);
testinds = cell(k,1);

train_size = ceil(opts.cvsplit*M);

for i=1:k
   traininds{i} = sort(inds{i}(1:train_size));
   testinds{i} = sort(inds{i}((train_size+1):end));
end

pool_train = parallel.pool.Constant(traininds);
pool_test = parallel.pool.Constant(testinds);
pool_A = parallel.pool.Constant(A);
pool_data = parallel.pool.Constant(data);
pool_Aeq = parallel.pool.Constant(Aeq);
pool_beq = parallel.pool.Constant(beq);
pool_B = parallel.pool.Constant(B);
pool_opts = parallel.pool.Constant(opts);


parfor i = 1:k
    wout = w0;
    
    tmp_scores = zeros(1, length(lambdas));
    
    for j = 1:length(lambdas)
        [ll, wout, exitflag] = cv(i, lambdas(j), wout, pool_train, ...
            pool_test, pool_A, pool_data, pool_Aeq, pool_beq, pool_B, pool_opts);
        
        if exitflag ~= 1 && exitflag ~= 2
            send(data_queue2, sprintf('Optimization may not have converged for lambda = %f (flag %d)', lambdas(j), exitflag));
        end
        
        tmp_scores(j) = ll;
        
        send(data_queue, 3);
    end
    
    scores(i, :) = tmp_scores;
end

[~, minlambdaind] = min(mean(scores, 1));
minlambda = lambdas(minlambdaind);

w = descent(w0, A, data.I, data.sigma, data.c, data.b, Aeq, beq, B, minlambda, opts.optimality, opts.steptol);
mapsi_waitbar(3)

% ------------------------- do final optimization ----------------------
minlambdaindse = find( mean(scores, 1) < (mean(scores(:,minlambdaind)) + std(scores(:,minlambdaind))/sqrt(size(scores, 1))), 1 );
minlambdase = lambdas(minlambdaindse);

wse = descent(w0, A, data.I, data.sigma, data.c, data.b, Aeq, beq, B, minlambdase, opts.optimality, opts.steptol);
mapsi_waitbar(3)

mapsi_waitbar(4)

% restore warnings
warning(orig_warnings);

% ------------------------- do error estimation -----------------------

[ wcovse ] = mapsi_mcmc( A, wse, minlambdase, data, opts);

end

% ======================== helper functions =============================

function [ll, wout, exitflag] = cv(i, l, w, pool_train, pool_test, ...
        pool_A, pool_data, pool_Aeq, pool_beq, pool_B, pool_opts)
% cross validation helper function
train = pool_train.Value{i};
test = pool_test.Value{i};

I = pool_data.Value.I;
sigma = pool_data.Value.sigma;
c = pool_data.Value.c;
b = pool_data.Value.b;

[wout, exitflag] = descent(w, pool_A.Value(train,:), I(train,:), ...
    sigma(train,:), c(train,:), b(train,:), pool_Aeq.Value, pool_beq.Value, pool_B.Value, ...
    l, pool_opts.Value.optimality, pool_opts.Value.steptol);

tmp = (I(test,:) - c(test,:) .* (pool_A.Value(test,:)*wout) - b(test,:)) ./ sigma(test,:);
ll = 0.5*mean(tmp.^2, 1);
end

% performs the optimization
function [w, exitflag] = descent( w0, A, I, sigma, c, b, Aeq, beq, B, lambda, opttol, steptol )

H = ( A .* ( c ./ sigma ) .^ 2 )' * A / length( I ) + lambda * B;

% ensure H is symmetric so quadprog doesn't complain
H = (H+H')/2;

f = mean(((b-I).*c./sigma.^2) .* A, 1)';

lb = zeros(size(w0));

opts = optimoptions('quadprog',...
    'Algorithm', 'interior-point-convex', 'Display', 'off',...
    'OptimalityTolerance',opttol,'StepTolerance',steptol,'MaxIterations',10000);

[w,~,exitflag,~] = quadprog(H,f,[],[],Aeq,beq,lb,[],w0,opts);

end

function inds = crossvalinds(N, k)
% returns a kx1 cell array of random permuatations of indices

inds = cell(k, 1);

for i=1:k
    inds{i} = randperm(N)';
end

end

function arg = numworkers()
p = gcp('nocreate');
if isempty(p)
  arg = 0;
else
  arg = p.NumWorkers;
end
end