function [quadratic_weights,thresholded_vicinity,linear_weights] = calculate_distanceweights(sub_xyz,group_xyz,thresh_mm)
% Calculates for all points on the group grid the normalized weights based on distance
% @param sub_xyz is the individuals mapping grid
% @param group_xyz the mapping grid for the whole group
% @param thresh_mm is the distance threshold: default threshold is 50 mm
% @returns quadratic_weigths as matrix of normalized weights (squared distance) 
% @returns utresholded_weigths as matrix of normalized weights (linear distance) 
% @returns thresholded_vicinity as matrix of grid points in vicinity of

if nargin <3, thresh_mm = 50; end  

%calculate distance for all points in the group-grid
all_dist                = pdist2(sub_xyz,group_xyz)';

if nargout>=1,
    t_quad_invdist          = (1./(all_dist.^2));
end

if nargout>=3,
    t_linear_invdist        = (1./all_dist);
end


% Calculate quadratic Weights based on inverse of quadratic distance
if nargout>=1,
    t_norm                  = repmat(sum(t_quad_invdist,2),1,size(t_quad_invdist,2));
    quadratic_weights       = (t_quad_invdist./t_norm);
    clear t_norm t_quad_invdist
end

% Calculate thresholded Weights based on distance for Latency 
% thresholded for maximal accepted distance
if nargout>=2,
    thresholded_vicinity    = all_dist<=thresh_mm;
end

% Calculate unthresholded linear Weights based on inverse of distance
% Consider that weights should be a weighted average, i.e. norm(x,1) should be equal to 1;
if nargout>=3,
    t_norm                  = repmat(sum(t_linear_invdist,2),1,size(t_linear_invdist,2));
    linear_weights          = (t_linear_invdist./t_norm);
    clear t_norm t_linear_invdist
end
