function CoG = calculate_CoG(mapping)

    base_val                = mapping.amp;    
    weights                 = (base_val./sum(base_val));
    CoG                     = weights'*mapping.xyz;
    CoG                     = single(CoG);
    
end