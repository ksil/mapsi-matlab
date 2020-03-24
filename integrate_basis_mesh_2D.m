function A = integrate_basis_mesh_2D( kernel, simp, in_verts, out_nodes, n_quad, status )

% kernel_func( x, y ) - a function that takes as
% input the data input domain x and output domain y.  It stores the angles phi and theta in 
% the respective columns of x.  kernel will be integrated over each local basis
% function on the surface of the sphere

% out_nodes - points at which the kernel is evaluated
% in_verts - list of coordinates for the vertices of different grid nodes
% simp - list of indices for the vertices connected to each simplex
% simplices must tile the plane cos_theta = [ -1, 1 ], phi = [ 0, 2 * pi ]
% n_quad - the number of quadrature points to use for numerical integration
% over the simplices

% data queue for logging progress
data_queue = parallel.pool.DataQueue;
afterEach(data_queue, @mapsi_waitbar);

if nargin < 6
    status = 0;
end

sz = size( out_nodes, 1 );
A = zeros( sz, size( in_verts, 1 ) );

memory_footprint_gb = 5 * sz * n_quad ^ 2 / 10 ^ 12 * 8; % 10^12 elements is about 8 GB

n_batches = ceil( memory_footprint_gb / 4 ); % Each batch should use roughly 4 GB of memory

if sz > numworkers()
    n_batches = max([n_batches, numworkers()]);
end

if status == 1
    mapsi_waitbar(2, size(simp,1)*n_batches);
end

batch_sz = floor(sz/n_batches);

A_cells = cell(n_batches,1);

parfor batch = 1:n_batches
    
    batch_start = batch_sz*(batch-1) + 1;
    batch_end = batch_sz*batch;
    
    if batch == n_batches
        batch_end = sz;
    end
    
    batch_indices = batch_start:batch_end;
    this_batch_sz = length(batch_indices);
    
    tmp_A = zeros(this_batch_sz, size(in_verts,1));
    
    % Expand the wave vectors to match dimensions of quadrature
    active_out_nodes = repelem( out_nodes( batch_indices, : ), n_quad ^ 2, 1 );
    
    
    % @x1, z = 1; @x2, z = 0; @x3, z = 0; returns z for x, y (vectorized)
    plane_func = @( x, y, x1, x2, x3 ) ( y * ( x2( 1 ) - x3( 1 ) ) + x * ( x3( 2 ) - x2( 2 ) ) + x3( 1 ) * x2( 2 ) - x2( 1 ) * x3( 2 ) ) ...
        / ( ( x2( 1 ) - x3( 1 ) ) * x1( 2 ) + x1( 1 ) * ( x3( 2 ) - x2( 2 ) ) + x3( 1 ) * x2( 2 ) - x2( 1 ) * x3( 2 ) ); 


    for i = 1:length( simp )

        V = in_verts( simp( i, : ), : );

        % fix the poles
        ind = find( abs(V( :, 2 )) < 1e-8 );
        if ( ~isempty( ind ) )
            not_ind = setdiff( [ 1 2 3 ], ind );
            V( ind, 1 ) = V( not_ind( 1 ), 1 );
        end;

        ind = find( abs(V( :, 2 ) - pi) < 1e-8 );
        if ( ~isempty( ind ) )
            not_ind = setdiff( [ 1 2 3 ], ind );
            V( ind, 1 ) = V( not_ind( 1 ), 1 ); 
        end;    

        % fix the prime meridian
        ind = find( abs(V( :, 1 )) < 1e-8 );
        if ( ~isempty( ind ) )
            if( max( V( :, 1 ) > pi ) )        
                V( ind, 1 ) = 2 * pi;
            end;
        end;    

        % Get quadrature rule on the simplex
        [ X, Y, Wx, Wy ] = triquad( n_quad, V );

        % Expand the quadrature for each output node
        X = repmat( X( : ), this_batch_sz, 1 );
        Y = repmat( Y( : ), this_batch_sz, 1 );   

        % Get the kernal part of the integrand
        f_save = kernel( [ X, Y ], active_out_nodes ) .* repmat( repmat( Wx, n_quad, 1 ) .* repelem( Wy, n_quad, 1 ), this_batch_sz, 1 ) .* sin( Y );

        % Loop over each vertex to get coefficient matrix contributions
        for j = 1:3   

            a = j; % Vertex with value 1 on the plane
            b = mod( j, 3 ) + 1; % Vertex with value 0 on the plane
            c = mod( j + 1, 3 ) + 1; % Vertex with value 0 on the plane    

            % Calculate the integrand as the product with the appropriate plane
            fA = f_save .* plane_func( X, Y, V( a, : ), V( b, : ), V( c, : ) );
            % Reshape the matrix and then add to the appropriate column of A
            fA = reshape( fA, n_quad ^ 2, this_batch_sz )';
            tmp_A( : , simp( i, j ) ) = tmp_A( : , simp( i, j ) ) + sum( fA, 2 );

        end;
        
        if status == 1
            send(data_queue, 3);
        end

    end;
    
    A_cells{batch} = tmp_A;
    
end;

% combine cell array from parallel for loop
cur = 0;
for i = 1:n_batches
    tmp_A = A_cells{i};
    A(cur + (1:size(tmp_A,1)),:) = tmp_A;
    cur = cur + size(tmp_A,1);
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