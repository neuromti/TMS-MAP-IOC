function IsArtifacted = find_Artifacted_IOC(ioc)
% artifacts (indicated by negative values)       
    
exclude_flag                = false(size(ioc.lat,1),1);
IsArtifacted                = (ioc.amp<0) | (ioc.lat <0) | exclude_flag;

end

