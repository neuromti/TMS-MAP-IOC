function xyz = get_DesignGridOrigin()    
    % Based on subject average:
    % load([folder.results.stats,'map_subject_data.mat'],'SUB')
    % GOxyz = []; GO = cat(3,SUB.GridOrigin); for xyz = 1:3, GOxyz = [GOxyz,nanmean(grpstats(squeeze(GO(:,xyz,:)),SUBID))]; end
    % xyz = GOxyz;
    xyz =  [-37.1589   -7.6858   68.9553];
end