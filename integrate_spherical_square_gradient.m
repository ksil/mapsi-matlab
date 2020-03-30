function B = integrate_spherical_square_gradient(verts, simp)
% calculates the matrix D such that w'*D*w yields the integral of the
% squared surface gradient over the surface of the sphere.
% Generally, a value of N leads to better accuracy for smooth functions.
%
% INPUT
% verts - vertices of mesh (Nx3)
% simp - list of simplices
% N - number of subdivisions for the trapezoidal method
%
% OUTPUT
% B - matrix of size Npoints x Npoints

Np = size(verts, 1);
B = zeros(Np, Np);

parfor t = 1:size(simp,1)
    V = verts(simp(t,:), :);            % vertices of each simplex
        
    p1 = simp(t,1);
    p2 = simp(t,2);
    p3 = simp(t,3);
        
    % calculate 3x3 contribution from these vertices
    Blocal = integrate_triangle(V);
    
    % create sparse matrix to efficiently add contribution to total D
    I = repmat([p1; p2; p3], 3, 1);
    J = repelem([p1; p2; p3], 3);
    Bsparse = sparse(I, J, Blocal(:), Np, Np);
    
    B = B + Bsparse;
end

end

function val = integrate_triangle(V)
% recursively divide the spherical triangle into other triangles
% doesn't touch the matrices V or W until the base case (should be pass by
% reference essentially, which is good)
%
% INPUT
% V - 3x3 matrix; each row is the (x,y,z) indices of a mesh vertex
%
% OUTPUT
% D - 3x3 local matrix used to calculate square gradient on this triangle

v1 = V(1,:) - V(3,:);
v2 = V(2,:) - V(3,:);

nv1 = norm(v1);
v1 = v1 / nv1;          % normalize v1
d12 = dot(v1,v2);

% estiamate surface gradient on triangle
Vinv = inv([nv1 0; d12 norm(v2 - d12*v1)]);

% need to multiply by B to acount for triangles shrinking
Y = [Vinv -sum(Vinv, 2)];

val = spherical_triangle_area(V(1,:), V(2,:), V(3,:)) * (Y'*Y);

end