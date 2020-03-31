function plot_dist(verts, simp, w, ec, plot_sub)
% PLOT_DIST plots the resulting MAPSI distribution on the surface of the
% sphere
%
% INPUT
% verts - vertices of mesh (in 3D cartesian coordinates)
% simp - list of simplices that comprise mesh
% w - weights of distribution at each vertex
% [ec] - edgecolor of mesh (default is none)
% [plot_sub] - if triangles should be subdivided and interpolated to
%   provide a smoother looking plot (default is 0)

if nargin < 4
    ec = 'none';
end

if nargin < 5
    plot_sub = 0;
end

if plot_sub > 0
    [verts, simp, w] = subdivide_interpolate(verts, simp, plot_sub, w);
end

figure
trisurf(simp, verts(:,1), verts(:,2), verts(:,3), w, 'EdgeColor', ec, 'FaceColor', 'interp');
axis square;
colorbar;
colormap(jet(256));

end