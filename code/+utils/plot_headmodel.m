function plot_headmodel(headmodel,V,CLIM)

if nargin<3
    CLIM = 2; %abs of upper and lower limit
    
end
FaceColor   = vec2headmodel(headmodel,V,CLIM); 
        
figure
set(gcf,'Position',[100,100,800,600])
ft_plot_mesh(headmodel,'edgealpha',0,'facecolor',FaceColor)
viewpoint   = [-90,75];
view(viewpoint)
[vx,vy,vz] = sph2cart(-90,90,128);
lightpos    =  [vx,vy,vz];
zoom(1.40)
light('Position',lightpos,'style','infinite');
material dull       	
lighting gouraud

end

function [FaceColor,TriValue] = vec2headmodel(headmodel,GridValue,CLIM)
   
Gridxyz  = utils.get_DesignGrid();

[quad_weights,threshed_weights] = utils.calculate_distanceweights(Gridxyz,headmodel.pos);
Vqery                           = utils.perform_weighting(GridValue,quad_weights,threshed_weights);   
Vqery                           = Vqery(:,2); %take only thresholded values

TriValue                        = (nanmean(Vqery(headmodel.tri),2));
BoundValue                      = BoundaryLimit(TriValue,CLIM);   
FaceColor                       = val2col(BoundValue);
    
    
function BoundValue = BoundaryLimit(TriValue,CLIM)

    absTriValue                     = abs(TriValue);
    sgnTriValue                     = sign(TriValue);    
    absTriValue(absTriValue>CLIM)   = CLIM;    
    BoundValue                      = (absTriValue.*sgnTriValue)./CLIM;
end
   
    
function FaceColor = val2col(BoundValue)   
    
    CMAP                        = get(groot,'DefaultFigureColormap');
    C_BIAS                      = (size(CMAP,1)./2-1);
    ColorIdx                    = 1+int16((C_BIAS+(BoundValue.*C_BIAS)));

    FaceColor                   = CMAP(ColorIdx,:);

end



end