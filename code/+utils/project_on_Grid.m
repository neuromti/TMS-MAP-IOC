% Core Function
% @param mapping with fields amp, lat, mep and xyz 
% @returns Vout with weighted/thresholded values projected from the actually employed to the predefined Mapping-Grid
function Vout = project_on_Grid(mapping)

Grid = get_GridPoints(mapping.xyz);
[quad_weights,threshed_weights] =  utils.calculate_distanceweights(Grid.raw.xyz,Grid.ave.xyz);
    
tmp.amp     = fix_on_Grid(mapping,'amp');
tmp.lat     = fix_on_Grid(mapping,'lat');
tmp.mep     = fix_on_Grid(mapping,'mep');
tmp.xyz     = Grid.raw.xyz;
Vout        = utils.perform_weighting(tmp,quad_weights,threshed_weights);

end

%% Supporting Functions----------------------------------------------------

% Function to derive Average Grid-Points from measurements for later
% projection on proper rectangular grid
% @params xyz as Matrix with xyz as column entries
% @returns Grid with fields for units and dimensions
% @returns Grid with fields for measured and average-proper xyz coordinates as Matrix with xyz as column entries
function Grid = get_GridPoints(xyz)

% Study Notebook: Maps of the left M1 and NPMA were acquired using a 7x15-point grid (5mm between points),
% centered 1cm anterior of M1-Hotspot.  Each point was stimulated 3 times with TMS.
Grid            = struct();
Grid.unit       = 'mm';

% Set Grid Parameters
Grid.X_Size     = int8(7); 
Grid.Y_Size     = int8(15);
Grid.Rep_Size   = int8(3);  

% Define Grid based on actually used grid-points
rawX                = (squeeze(nanmean(reshape(xyz(:,1),Grid.Rep_Size,Grid.Y_Size,Grid.X_Size),1)));
rawY                = flipEven(squeeze(nanmean(reshape(xyz(:,2),Grid.Rep_Size,Grid.Y_Size,Grid.X_Size),1)));
rawZ                = flipEven(squeeze(nanmean(reshape(xyz(:,3),Grid.Rep_Size,Grid.Y_Size,Grid.X_Size),1)));
Grid.raw.xyz        = demesh(rawX,rawY,rawZ);

aveX                = repmat(nanmean(rawX,1),Grid.Y_Size,1);
aveY                = repmat(nanmean(rawY,2),1,Grid.X_Size);
aveZ                = repmat(nanmean(nanmean(rawZ)),Grid.Y_Size,Grid.X_Size);
Grid.ave.xyz        = demesh(aveX,aveY,aveZ);

end

% Function to de-mesh XYZ
% @params   3 meshgrids for X Y Z
% @returns  Matrix with xyz as column entries
function xyz = demesh(x,y,z)
    xyz = cat(1,reshape(x,1,[]),reshape(y,1,[]),reshape(z,1,[]))';
end

% Control Function
% - Averages the repeated measures for each grid point 
% - Fixes NaNs by nearest-neighbour interpolation
% @params mapping with fields amp, lat, mep and xyz 
% @params field with string of fieldname to perform function on
% @returns Vout, the avreaged and repaired values ordered according to xyz
function Vout = fix_on_Grid(mapping,field)

eval(['rawV = mapping.',field,';']);
V       = squeeze(nanmean(reshape(rawV,3,15,7),1));
V       = fix_NaN(V);
Vout    = reshape(flipEven(V),1,[])';

end

% Function to fix NaNs by average of neighbours 
% @params Vin as meshgrid of values
% @returns Vout as meshgrid of values
function Vout = fix_NaN(Vin)

Vout        = Vin;
[X,Y]       = meshgrid(1:size(Vin,2),1:size(Vin,1));
[Xq,Yq]     = meshgrid(0:size(Vin,2)+1,0:size(Vin,1)+1);
Vq          = interp2(X,Y,Vin,Xq,Yq,'linear');

for i_x = 1:size(Vin,1)
    for i_y = 1:size(Vin,2)
        if isnan(Vin(i_x,i_y))
            Vout(i_x,i_y) = nanmean(nanmean((Vq(i_x:i_x+2,i_y:i_y+2))));            
        end
    end
end
end

% Function to flip even meshgrid columns for proper alignment 
% Introduced as mapping-grid was recorded oscillatory (i.e. up-down-up-down....)
% @params raw unflipped meshgrid
% @returns flipped Y-coordinates
function y = flipEven(y)

for ix = 1:size(y,2)
    if ~mod(ix,2) % if even
        y(:,ix)  = flipud(y(:,ix));
    end
end
end