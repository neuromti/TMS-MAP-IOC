function [grid_xyz,map_dimensions,X,Y,Z] = make_groupgrid(lcl_folder,IsArtifacted,D)
% Define the outlines of the group-wise grid for projection of individual
% maps
% @return grid_xyz, specifying the coordinates of the group-wise grid in (X,Y,Z)
% @return 
if nargin <3,
    D           = dir(lcl_folder);
    D([1,2])    = [];
end

% extension of the grid by 2 cm
grid_extension  = 2;

XYZ  = [];
for idx_dataset = 1:length(D),
    if ~IsArtifacted(idx_dataset),
        load([lcl_folder,D(idx_dataset).name]);    
        XYZ = cat(1,XYZ,cat(1,max(mapping.xyz),min(mapping.xyz)));
    end
end
 
xyz_upperbound  = fix(max(XYZ)+grid_extension);
xyz_lowerbound  = fix(min(XYZ)-grid_extension);
[X,Y,Z]         = meshgrid(xyz_lowerbound(1):xyz_upperbound(1),xyz_lowerbound(2):xyz_upperbound(2),xyz_lowerbound(3):xyz_upperbound(3));    
map_dimensions  = size(X);
grid_xyz        = cat(1,reshape(X,1,[]),reshape(Y,1,[]),reshape(Z,1,[]))';