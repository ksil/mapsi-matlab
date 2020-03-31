function plot_dist_error(verts, simp, wcov, ec, plot_sub)
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

diag_wcov = diag(wcov);

if plot_sub > 0
    [verts, simp, diag_wcov] = subdivide_interpolate(verts, simp, plot_sub, diag_wcov);
end

figure
cmap = [1 1 1;
        2/3 2/3 1;
        1/3 1/3 1;
        0 0 1;
        1/3 0 2/3;
        2/3 0 1/3;
        1 0 0];
    
trisurf(simp, verts(:,1), verts(:,2), verts(:,3), log10(abs(diag_wcov)), 'EdgeColor', ec,'FaceColor','interp');
axis square;
cbar = colorbar;
caxis([-8, -1])
colormap(cmap)
ylabel(cbar, 'log_{10}|Cov(w_i, w_j)|')

end