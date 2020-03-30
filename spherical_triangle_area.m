function val = spherical_triangle_area(v1, v2, v3)
    % returns the area of the triangle on a unit sphere with vertices
    % v1, v2, and v3
    
    n1 = cross(v1, v2); n1 = n1/norm(n1);
    n2 = cross(v2, v3); n2 = n2/norm(n2);
    n3 = cross(v3, v1); n3 = n3/norm(n3);
    
    % ensure outward facing normals
    if dot(n1, v3) > 0
        n1 = -n1; n2 = -n2; n3 = -n3;
    end
    
    % calculate inner angles
    a1 = pi - acos(max(min(dot(n1, n2), 1.0), -1.0));
    a2 = pi - acos(max(min(dot(n2, n3), 1.0), -1.0));
    a3 = pi - acos(max(min(dot(n3, n1), 1.0), -1.0));
    
    val = a1 + a2 + a3 - pi;
end