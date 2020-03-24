function [ wcov ] = mapsi_mcmc( A, w, minlambda, data, opts)
% MAPSI_MCMC performs adaptive Monte Carlo sampling to estimate the
% covariance matrix associated with the weights of the size distribution
%
% INPUT
% A - the mapsi kernel matrix
% w - the MAP estimate of the weights
% minlambda - the lambda chosen from cross validation
% data - a mapsi_data struct
% opts - a mapsi_options struct
%
% OUTPUT
% wcov - covariance matrix of the weights

orig_warnings = warning;
warning('off','all')

mapsi_waitbar(0);

if isempty(gcp('nocreate'))
    parpool;
    fprintf('\n')
end

trials = opts.mctrials;
rng('shuffle');

if isrow(w)
    w = w';
end

[in_verts, simp] = create_kurihara_mesh(opts.N);
Aeq = integrate_basis_mesh_2D(@(in,out) 1, simp, in_verts, 0, opts.n_quad);
B = integrate_square_gradient( in_verts, simp, opts.n_quad );

% define matrices for objective function
H = ( A .* ( data.c ./ data.sigma ) .^ 2 )' * A + minlambda * length(data.I) * B;
f = sum(((data.b-data.I).*data.c./data.sigma.^2) .* A, 1);
H = (1./Aeq') .* H .* (1./Aeq);
f = f ./ Aeq;

H = parallel.pool.Constant(H);
f = parallel.pool.Constant(f);

% adaptively determine best value of kappa via pseudo-binary search
% aim for acceptance rate between 0.35 and 0.4
kappa = 0;
kappas = [1 0] / length(w);
accepted = [0 1];
found = false;

for i=1:500
    eval = parfeval(@mc, 3, H, f, Aeq, w, 1e3, kappas(1), 0);
    [~, ~, ~, accepted(1)] = fetchNext(eval);
    
    if accepted(1) < 0.35
        kappas(1) = mean(kappas); % double
    elseif accepted(1) > 0.4
        kappas(2) = kappas(1);
        accepted(2) = accepted(1);
        kappas(1) = kappas(1) * 2;
    else
        kappa = kappas(1);
        found = true;
        break
    end
end

if ~found
    error('Could not achieve desired acceptance rate')
end

mapsi_waitbar(1, 'Estimating error')

workers = numworkers();
trialsperworker = ceil(trials/workers);
tottrials = trialsperworker*workers;

samples = zeros(length(w), tottrials);
lls = zeros(1, tottrials);
totaccepted = 0;

for i=1:workers
    evals(i) = parfeval(@mc, 3, H, f, Aeq, w, ceil(trials/workers), kappa, 1);
end

mapsi_waitbar(2, tottrials);

trialscompleted = 0;
while trialscompleted < tottrials
    trialscompleted = 0;
    for i = 1:workers
        trialscompleted = trialscompleted + ceil(trialsperworker/100)*length(evals(i).Diary);
    end
    
    if trialscompleted > tottrials
        trialscompleted = tottrials;
    end
    
    mapsi_waitbar(3, trialscompleted);
    
    pause(0.05);
end

for i=1:workers
    [i, s, l, a] = fetchNext(evals);
    totaccepted = totaccepted + a/workers;
    samples(:, (1 + (i-1)*trialsperworker):(i*trialsperworker)) = s;
    lls((1 + (i-1)*trialsperworker):(i*trialsperworker)) = l;
end

mapsi_waitbar(1, sprintf('MC Acceptance: %.3f', totaccepted));

% calculate covariance matrix
wcov = cov(samples');

mapsi_waitbar(4)

% restore warnings
warning(orig_warnings);
end

function neww = sample_w(w, kappa, dims)
v = sqrt(w);

% sample in n dimensional hypersphere
neww = randn(dims, 1);

% project onto tangent space
neww = neww - v*dot(v,neww);
neww = neww * kappa * (rand)^(1/(dims-1)) / norm(neww);

% stereographically project onto original hypersphere
neww = v + neww;
neww = neww .* neww;
neww = neww / sum(neww);

end

function [samples, lls, accepted] = mc(H, f, Aeq, w0, trials, kappa, shouldprint)
H = H.Value;
f = f.Value;
cur = Aeq' .* w0;
dims = length(w0);
samples = zeros(dims, trials);
lls = zeros(1, trials);

llcur = -(cur'*H*cur)/2 - f*cur;

i = 1;
accepted = 0;
while i <= trials
    next = sample_w(cur, kappa, dims);
    
    llnext = -(next'*H*next)/2 - f*next;
    
    % record sample
    if (llnext - llcur) > log(rand)
        samples(:, i) = next ./ Aeq';
        cur = next;
        llcur = llnext;
        accepted = accepted + 1;
    else
        samples(:, i) = cur ./ Aeq';
    end
    
    lls(i) = llcur;
    
    i = i+1;
    
    if shouldprint && mod(i, ceil(trials/100)) == 0
        fprintf('.')
    end
end
if shouldprint
    fprintf('.')
end

accepted = accepted/trials;
end

function arg = numworkers()
p = gcp('nocreate');
if isempty(p)
  arg = 0;
else
  arg = p.NumWorkers;
end
end