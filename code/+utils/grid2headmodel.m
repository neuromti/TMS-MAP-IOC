function FaceColor = grid2headmodel(headmodel,V)

   
design_xyz                      = utils.get_DesignGrid();
[quad_weights,threshed_weights] = utils.calculate_distanceweights(design_xyz,headmodel.pos);
Vqery                           = utils.perform_weighting(V,quad_weights,threshed_weights);   
Vqery                           = Vqery(:,2); %take only thresholded values
Tqery                           = (nanmean(Vqery(headmodel.tri),2));
 

CMAP                        = get(groot,'DefaultFigureColormap');
CLIM                        = 2;
C_BIAS                      = (size(CMAP,1)./2-1);
ColorIdx                    = abs(Tqery./CLIM);
ColorIdx(ColorIdx>1)        = 1;
ColorIdx                    = int16(C_BIAS+(ColorIdx.*C_BIAS));

FaceColor                   = CMAP(ColorIdx,:);
