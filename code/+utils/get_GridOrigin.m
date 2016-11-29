function grid_origin = get_GridOrigin(mapping)
    
    UseDefinedOrigin = true;

    if UseDefinedOrigin 
        % Study Notebook: Maps of the left M1 and NPMA were acquired using a 7x15-point grid (5mm between points),
        % centered 1cm anterior of M1-Hotspot.  Each point was stimulated 3 times with TMS.        
        xyz                 = squeeze(nanmean(reshape(mapping.xyz,3,15,7,3),1));
        grid_origin         = squeeze(xyz(8,4,:))';     
    
    else
        grid_origin         = mean(mapping.xyz(1:315,:));
    end
    
end