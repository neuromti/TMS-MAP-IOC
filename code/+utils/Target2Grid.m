function V = Target2Grid(XYZ)
    W = [];
    for run_idx = 1:size(XYZ,1),
        Grid        = double(utils.get_DesignGrid);
        w           = 1./(pdist2(Grid,XYZ(run_idx,:)).^2);
        w           = min(w,1);
        w           = w./sum(w);
        W           = cat(2,W,w);
    end
    V           = nanmean(W,2);
    V           = V./max(max(V));
end