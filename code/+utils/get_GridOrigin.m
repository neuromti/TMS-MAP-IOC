function grid_origin = get_GridOrigin(mapping)
    grid_origin  = mean(mapping.xyz(1:315,:));
end