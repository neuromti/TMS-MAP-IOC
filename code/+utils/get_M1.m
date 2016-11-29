function xyz = get_M1()
    % According to Mayka, M. A., Corcos, D. M., Leurgans, S. E., & Vaillancourt, D. E. (2006):
    % Three-dimensional locations and boundaries of motor and premotor cortices as defined by functional brain imaging: a meta-analysis. 
    % NeuroImage, 31(4), 1453–1474. https://doi.org/10.1016/j.neuroimage.2006.02.004
    % Consider that any MNI coordinates not reported in Talairach space were converted using the transformation equations for above the AC line (z > 0):
    % xV = 0.9900x yV = 0.9688y + 0.0460z, zV = -0.0485y + 0.9189z 
    % xyz_in_Tailarach    = [-37, -21, 58]; 
    % Tailarach2MNI       = [0.9900,0,0;0,0.9688,0.0460;0,-0.0485,0.9189];
    % xyz_in_MNI          = (Tailarach2MNI*xyz_in_Tailarach')'; 
    % Result as literals for speed
    xyz  =  [-36.6300  -17.6768   54.3147];
    
    
end

