function [xyz,name] = get_landmark(name,coordsys)

% According to Mayka, M. A., Corcos, D. M., Leurgans, S. E., & Vaillancourt, D. E. (2006):
% Three-dimensional locations and boundaries of motor and premotor cortices as defined by functional brain imaging: a meta-analysis. 
% NeuroImage, 31(4), 1453–1474. https://doi.org/10.1016/j.neuroimage.2006.02.004
% Consider that any MNI coordinates not reported in Talairach space were converted using the transformation equations for above the AC line (z > 0):
% xV = 0.9900x yV = 0.9688y + 0.0460z, zV = -0.0485y + 0.9189z 
% xyz_in_Tailarach    = [-37, -21, 58]; 
% Tailarach2MNI       = [0.9900,0,0;0,0.9688,0.0460;0,-0.0485,0.9189];
% xyz_in_MNI          = ([0.9900,0,0;0,0.9688,0.0460;0,-0.0485,0.9189]*xyz_in_Tailarach')'; 

tailarach2mni      = @(x)([0.9900,0,0;0,0.9688,0.0460;0,-0.0485,0.9189]*x')';
%
%MAYKA.MPMC        = [-2 -1 54];
MAYKA.preSMA      = [-3 6 53];
MAYKA.SMAproper   = [-2 -7 55];
%MAYKA.LPMC        = [-26 -6 56]; 
MAYKA.PMd         = [-30 -4 58];
MAYKA.PMv         = [-50 5 22];
%MAYKA.SMC         = [-39 -21 54];
MAYKA.M1          = [-37 -21 58];
MAYKA.S1          = [-40 -24 50];

CoordSys_enum       = {'MNI','TAILARACH'};
if nargin <2, coordsys = 'MNI'; end
if nargin <1, name = fieldnames(MAYKA); end
if ~any(ismember(fieldnames(MAYKA),name)), error('LANDMRK:Unspecified','Landmark not correctly specified'); end %case sensitive

xyz = [];
for name_idx = 1 : length(name)
    
if strcmpi(CoordSys_enum{1},coordsys) %MNI
    xyz = cat(1,xyz,tailarach2mni(getfield(MAYKA,name{name_idx})));
elseif strcmpi(CoordSys_enum{2},coordsys) %TAILARACH
    xyz = cat(1,xyz,getfield(MAYKA,name{name_idx}));
else
    error('LANDMRK:CoordSysUnspecified','Coordinate System not correctly specified'); 
end

end

end