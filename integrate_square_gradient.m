function B = integrate_square_gradient( in_verts, simp, n_quad  )

% integrate_square_gradient

% INPUTS:
% in_verts - list of coordinates for the vertices of different grid nodes
% simp - list of indices for the vertices connected to each simplex
% simplices must tile the plane cos_theta = [ -1, 1 ], phi = [ 0, 2 * pi ]
% n_quad - the number of quadrature points to use for numerical integration
% over the simplices

B = zeros( length( in_verts ), length( in_verts ) );

% @x1, z = 1; @x2, z = 0; @x3, z = 0; returns gradient of the plane
plane_gradient = @( x1, x2, x3 ) [ x3( 2 ) - x2( 2 ), x2( 1 ) - x3( 1 )] ...
    ./ ( ( x2( 1 ) - x3( 1 ) ) * x1( 2 ) + x1( 1 ) * ( x3( 2 ) - x2( 2 ) ) + x3( 1 ) * x2( 2 ) - x2( 1 ) * x3( 2 ) ); 

for i = 1:length( simp )
    
    vert_inds = simp( i, : );
    V = in_verts( vert_inds, : );     % all the vertices for this simplex

    % Fix the poles
    ind = find( V( :, 2 ) == 0 );
    if ( ~isempty( ind ) )
        not_ind = setdiff( [ 1 2 3 ], ind );
        V( ind, 1 ) = V( not_ind( 1 ), 1 );
    end;

    ind = find( V( :, 2 ) == pi );
    if ( ~isempty( ind ) )
        not_ind = setdiff( [ 1 2 3 ], ind );
        V( ind, 1 ) = V( not_ind( 1 ), 1 ); 
    end;    
    
    % Fix the prime meridian
    ind = find( V( :, 1 ) == 0 );
    if ( ~isempty( ind ) )
        if( max( V( :, 1 ) > pi ) )        
            V( ind, 1 ) = 2 * pi;
        end;
    end;
    
    
    % Get quadrature rule on the simplex
    [ X, Y, Wx, Wy ] = triquad( n_quad, V );   
    

    % Loop over each vertex to get coefficient matrix contributions
    for j = 1:3   

        a = j; % Vertex with value 1 on the plane
        b = mod( j, 3 ) + 1; % Vertex with value 0 on the plane
        c = mod( j + 1, 3 ) + 1; % Vertex with value 0 on the plane           
        
        grad_j = plane_gradient( V( a, : ), V( b, : ), V( c, : ) );
        
        for k = j:3
            
            if ( j == k )
            
                grad_k = grad_j;
                
            else
                
                a = k; % Vertex with value 1 on the plane
                b = mod( k, 3 ) + 1 ; % Vertex with value 0 on the plane
                c = mod( k + 1, 3 ) + 1; % Vertex with value 0 on the plane         
        
                grad_k = plane_gradient( V( a, : ), V( b, : ), V( c, : ) );
                
            end;
            
            % Calculate the integrand as the product with the appropriate plane
            % and weight by the appropriate phi weights
            fB = Wx' * ( ( grad_j( 1 ) * grad_k( 1 ) ./ sin( Y ) .^ 2 + grad_j( 2 ) * grad_k( 2 ) ) .* sin( Y ) ) * Wy;
           
            B( vert_inds( j ), vert_inds( k ) ) = B( vert_inds( j ), vert_inds( k ) ) + fB;
            
            if ( k > j )
        
                B( vert_inds( k ), vert_inds( j ) ) = B( vert_inds( k ), vert_inds( j ) ) + fB;
                
            end;
            
        end;
        
    end;
    
end;