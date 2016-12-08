function grid_origin = get_GridOrigin(mapping,UseDefinedOrigin)
    
if nargin <2,
    UseDefinedOrigin = true;
end

if UseDefinedOrigin 
    % Study Notebook: Maps of the left M1 and NPMA were acquired using a 7x15-point grid (5mm between points),
    % centered 1cm anterior of M1-Hotspot.  Each point was stimulated 3 times with TMS.   
    % Use Coordinates of this point for grid_origin
    xyz                 = squeeze(nanmean(reshape(mapping.xyz,3,15,7,3),1));
    grid_origin         = squeeze(xyz(8,4,:))';     

else % Use Coordinates of the average grid positions for grid_origin
    grid_origin         = mean(mapping.xyz(1:315,:));
end
    
end