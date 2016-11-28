function remove_flag = find_artifactedTrials(mapping)

exclude_flag                = true(size(mapping.xyz,1),1);
exclude_flag(1:315)         = false;
remove_flag                 = mapping.amp<0 | mapping.lat <0  | exclude_flag;

end

