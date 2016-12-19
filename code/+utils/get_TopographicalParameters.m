function [CoG,Spread,Area] = get_TopographicalParameters(member)

    Area    = sum(sum(member));
    xyz     = utils.get_DesignGrid(utils.mesh2vec(member));
    CoG     = mean(xyz,1);    
    Spread  = sqrt(mean(pdist2(xyz,CoG).^2));
  %
end