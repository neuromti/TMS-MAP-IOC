function xyz = get_DesignGridOrigin()    
    % Based on subject average:
    % load([folder.results.stats,'map_subject_data.mat'],'SUB')
    % GOxyz = []; GO = cat(3,SUB.GridOrigin); for xyz = 1:3, GOxyz = [GOxyz,nanmean(grpstats(squeeze(GO(:,xyz,:)),cat(1,SUB.subID)))]; end
    % xyz = GOxyz;
    xyz =  [-36.8820   -8.5896   69.3859];
end