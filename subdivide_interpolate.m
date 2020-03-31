function [verts, simp, w] = subdivide_interpolate(verts, simp, N, w)
% subdivide a mesh and interpolate the weights, w, at each point
%
% INPUT
% verts - vertices of mesh
% simp - list of simplices of mesh
% N - number of subdivisions
% [w] - weights of function at each mesh point
%
% OUTPUT
% subdivided and interpolated quantites of above

if N < 1    % don't need to subdivide then!
    return
end

if nargin < 4
    w = zeros(size(verts, 1), 1);
end

pair = @(k1,k2) max(k1,k2)*(max(k1,k2)-1)/2 + min(k1,k2);

% ------------ construct edges -------------------
nbors = zeros(size(verts, 1), 10);
num_nbors = zeros(size(verts, 1), 1);
for t = 1:size(simp, 1)
    v1 = simp(t,1); v2 = simp(t,2); v3 = simp(t,3);
    
    vpairs = [v1 v2 v3 v1];
    
    for i = 1:3
        % edge v1-v2
        m1 = min(vpairs(i), vpairs(i+1));
        m2 = max(vpairs(i), vpairs(i+1));
        found = false;
        
        % see if edge has already been added to list
        for j = 1:num_nbors(m1)
            if nbors(m1, j) == m2
                found = true;
                break
            end
        end
        
        if ~found
            num_nbors(m1) = num_nbors(m1) + 1;
            nbors(m1, num_nbors(m1)) = m2;
        end
    end
    
    if any(num_nbors > size(nbors,2) - 2)
        new_nbors = zeros(size(nbors,1), 2*size(nbors,2));
        new_nbors(:, 1:size(nbors,2)) = nbors;
        nbors = new_nbors;
    end
end

edges = zeros(sum(num_nbors), 2);

count = 1;
for i = 1:size(nbors, 1)
    for j = 1:num_nbors(i)
        edges(count,:) = [i nbors(i, j)];
        count = count + 1;
    end
end

% --------- create four triangles for every bigger triangle -----------
for i = 1:N
    Nsimp = size(simp,1);
    Nverts = size(verts,1);
    Nedges = size(edges,1);
    
    newverts = [verts; zeros(Nedges, 3)];    % one new vertex per edge
    neww = [w; zeros(Nedges, 1)];
    newsimp = zeros(4*Nsimp, 3);             % 4x the # of triangles
    newedges = zeros(2*Nedges + 3*Nsimp, 2);
    
    pairmap = java.util.HashMap; % map each edge to new vertex
    
    for j = 1:Nedges
        v1 = edges(j,1); v2 = edges(j,2);
        
        % add new vertex
        newvert = (verts(v1,:) + verts(v2,:)) / 2;
        
        % add new edges
        newverts(Nverts+j,:) = newvert;
        newedges(2*j-1,:) = [edges(j,1), Nverts+j];
        newedges(2*j,:) = [Nverts+j, edges(j,2)];
        
        % interpolate weights
        neww(Nverts+j) = (w(v1) + w(v2)) / 2;
        
        % for each edge, point to the new vertex created
        key = pair(v1, v2);
        pairmap.put(key, Nverts+j);
    end
    
    for k = 1:Nsimp
        v1 = simp(k,1); v2 = simp(k,2); v3 = simp(k,3);
        v12 = pairmap.get(pair(v1, v2));
        v13 = pairmap.get(pair(v1, v3));
        v23 = pairmap.get(pair(v2, v3));
        
        % add new edges
        newedges(2*Nedges + 3*k-2,:) = [v12 v13];
        newedges(2*Nedges + 3*k-1,:) = [v12 v23];
        newedges(2*Nedges + 3*k,:) = [v13 v23];
        
        % construct new triangles
        newsimp(4*k-3,:) = [v1 v12 v13];
        newsimp(4*k-2,:) = [v2 v12 v23];
        newsimp(4*k-1,:) = [v3 v13 v23];
        newsimp(4*k,:) = [v12 v13 v23];
    end
    
    % update with new matrices
    verts = newverts;
    w = neww;
    simp = newsimp;
    edges = newedges;
    clear pairmap
end

% normalize vertices to lie on sphere (important to do this now so that
% weights are interpolated correctly)
verts = verts ./ sqrt(sum(verts.*verts, 2));

end