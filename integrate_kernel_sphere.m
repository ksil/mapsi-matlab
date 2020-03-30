function A = integrate_kernel_sphere(kernel, simp, in_verts, out_nodes, n_sub, status)
% Calculates the kernel matrix over basis functions defined on the sphere
% over a triangular mesh
%
% INPUT
% kernel( x, y ) - a function that takes as
%   input the data input domain x and output domain y.  It stores the angles phi and theta in 
%   the respective columns of x.  kernel will be integrated over each local basis
%   function on the surface of the sphere
% simp - list of indices for the vertices connected to each simplex
% in_verts - list of coordinates for the vertices of different grid nodes
%   (in 3 dimensional Cartesian coordinates)
% out_nodes - points at which the kernel is evaluated
% n_sub - number of subdivisions of spherical triangles for integration
% [status] - 1 if waitbar should be printed
%
% OUTPUT
% A - the kernel matrix

% data queue for logging progress
data_queue = parallel.pool.DataQueue;
afterEach(data_queue, @mapsi_waitbar);

if nargin < 6
    status = 0;
end

sz = size(out_nodes, 1);
A = zeros(sz, size(in_verts, 1));

n_batches = numworkers();
batch_sz = floor(sz/n_batches);

if status == 1
    mapsi_waitbar(2, size(simp,1)*n_batches);
end

A_cells = cell(n_batches,1);

parfor batch = 1:n_batches
    
    % define batch of out nodes
    batch_start = batch_sz*(batch-1) + 1;
    batch_end = batch_sz*batch;
    
    if batch == n_batches
        batch_end = sz;
    end
    
    batch_indices = batch_start:batch_end;
    this_batch_sz = length(batch_indices);
    
    tmp_A = zeros(this_batch_sz, size(in_verts,1));
    
    % Define function to integrate over spherical triangles
    % should return a this_batch_sz x 3 matrix of values
    this_out_nodes = out_nodes(batch_indices, :);
    f = @(p,t,w,b) kernel([p t], this_out_nodes) .* [b(1) b(2) b(3)];
    
    % loop over each simplex
    for t = 1:size(simp, 1)
        % vertices of simplex
        V = in_verts(simp(t,:), :);

        % contribution to A from these 3 vertices
        contrib = integrate_on_spherical_triangle(f, V, [0;0;0], n_sub);
        
        % add contribution to tmp_A
        tmp_A(:, simp(t,:)) = tmp_A(:, simp(t,:)) + contrib;
        
        if status == 1
            send(data_queue, 3);
        end
    end
    
    A_cells{batch} = tmp_A;
    
end

% combine cell array from parallel for loop
for i = 1:n_batches
    batch_start = batch_sz*(i-1) + 1;
    batch_end = batch_sz*i;
    
    if i == n_batches
        batch_end = sz;
    end
    
    batch_indices = batch_start:batch_end;

    A(batch_indices, :) = A_cells{i};
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