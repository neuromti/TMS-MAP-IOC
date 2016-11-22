function [projected_grid] = interpolate_latency(in_vicinity,sub_lat)    

%remove grid points with artifacts or no MEP from raster
sub_lat(sub_lat<=0)             = NaN;
%
%take_val                        = (repmat(sub_lat,1,size(in_vicinity,1)))';

% Find grid-points in vicinity, and calculate weights based on weighted histogram of latencies in vicinity
[row_latency,column_grid]       = find(in_vicinity');    
lat_inAnyVicinity               = sub_lat(row_latency);

uniq_lat                        = unique(lat_inAnyVicinity);
uniq_lat(isnan(uniq_lat))       = [];
uniq_grid                       = unique(column_grid);

histo_mat                       = crosstab(lat_inAnyVicinity,column_grid);
weight_mat                      = histo_mat./repmat(sum(histo_mat),size(histo_mat,1),1); 
L                               = weight_mat'*uniq_lat;

% Remove grid-points without valid latencies
uniq_grid(isnan(L))             = [];
valid_L                         = L;
valid_L(isnan(L))               = [];

% Construct output-grid and center latencies to remove offset
projected_grid                  = zeros(size(in_vicinity,1),1);
projected_grid(uniq_grid)       = valid_L-mean(valid_L);