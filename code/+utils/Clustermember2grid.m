function V = Clustermember2grid(ClusterMember)
     
Grid    = utils.get_DesignGrid;
sel     = utils.mesh2vec(ClusterMember);
xyz     = double(Grid(sel,:));
V       = utils.Target2Grid(xyz);

end