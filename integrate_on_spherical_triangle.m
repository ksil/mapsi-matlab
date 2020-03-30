function val = integrate_on_spherical_triangle(f, V, W, N)
% integrate a function on a spherical triangle using the trapezoidal method
%
% INPUT
% f - a function handle of the form f(phi, theta, w, b) where
%   phi - azimuthal angle [0 2*pi]
%   theta - polar angle [0 pi]
%   w - weight at center of triangle (from linear interpolation of weights
%       at vertices)
%   b - barycentric coordinate at center of triangle
% V - 3x3 matrix; each row is the (x,y,z) indices of a mesh vertex
% W - 3 column vector of the weights at each mesh vertex
% N - number of subdivisions for the trapezoidal method
% s - the partial sum (should be 0 at the beginning

if ~iscolumn(W)
    W = W';
end

val = recursive_integrate(f, V, W, N, eye(3));

end

function val = recursive_integrate(f, V, W, N, B)
% recursively divide the spherical triangle into other triangles
% doesn't touch the matrices V or W until the base case (should be pass by
% reference essentially, which is good)
%
% INPUT
% same as above except
% B - 3x3 matrix of barycentric coordinates where each row is barycentric
%   coordinates of each vertex

if N == 0
    newV = B*V;      % doesn't matter if normalized;
    newW = B*W;
    
    [phi, theta, ~] = cart2sph(mean(newV(:,1)), mean(newV(:,2)), mean(newV(:,3)));
    theta = pi/2 - theta; % convert elevation to polar angle
    
    f_cent = f(phi, theta, mean(newW), mean(B, 1));
    
    val = f_cent * spherical_triangle_area(newV(1,:), newV(2,:), newV(3,:));
    return
end

% if 1 is bottom left and rotate counter clockwise:
% bottom left triangle
val =       recursive_integrate(f, V, W, N-1, [ B(1,:);
                                               (B(1,:) + B(2,:))/2;
                                               (B(1,:) + B(3,:))/2]);
% bottom right triangle
val = val + recursive_integrate(f, V, W, N-1, [(B(1,:) + B(2,:))/2;
                                                B(2,:);
                                               (B(2,:) + B(3,:))/2]);
% top triangle
val = val + recursive_integrate(f, V, W, N-1, [(B(1,:) + B(3,:))/2
                                               (B(2,:) + B(3,:))/2
                                                B(3,:)]);
% middle triangle
val = val + recursive_integrate(f, V, W, N-1, [(B(1,:) + B(3,:))/2
                                               (B(2,:) + B(3,:))/2;
                                               (B(1,:) + B(2,:))/2]);
end