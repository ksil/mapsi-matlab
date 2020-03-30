function [verts, simp] = kurihara_mesh(N)
% Create a Kurihara mesh with roughly equal surface area per simplex.
% The Kurihara mesh uses 4(N+1) equitorial vertices and 
%
% INPUTS:
% N - the number of equitorial nodes in one quadrant of the unit sphere.
%
% OUTPUTS:
% verts - vertices of the mesh (in 3 dimensions)
% simp - lists of indices for the vertices that comprise each simplex

% Look over lattitudes starting at the equator
k = 0;
for i = (N + 1):-1:0
    
    if i ~= 0
    
        % Place the right number of nodes around the pole at this lattitude
        for j = 1:( 4 * i )

            k = k + 1;
            phi(k) = (j - 1) / (4*i) * 2*pi;
            theta(k) = i / (N + 1) * pi / 2;

        end
    
    else
    
        % Place the pole
        k = k + 1;
        phi(k) = 0;
        theta(k) = 0;
        
    end
        
end

% Set to correct coordinates for MATLABs sph2cart
theta = pi/2 - theta;

% Duplicate the grid for the southern hemisphere
phi2 = phi((4*(N + 1) + 1):end);
theta2 = -theta((4*(N + 1) + 1):end);

[x, y, z] = sph2cart(phi, theta, ones(size(theta)));
[x2, y2, z2] = sph2cart(phi2, theta2, ones(size(theta2)));

verts = [x, x2; y, y2; z, z2]';

% Find the simplices by computing the convex hull of the Cartesian
% coordinates and compute the vertices in the phi [0,2pi], theta [0,pi]
% coordinate system.
simp = convhulln(verts);

end