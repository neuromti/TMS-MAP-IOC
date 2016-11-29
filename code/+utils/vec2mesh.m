function xyz = vec2mesh(xyz)
    if length(size(xyz))==2,
        xyz = reshape(xyz,15,7);
    elseif length(size(xyz))==3,
        xyz = reshape(xyz,15,7,size(xyz,2));
    else
        error('Function can only transform up to 3D-Matrices');
    end
        
end
