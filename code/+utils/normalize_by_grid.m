function  sub = normalize_by_grid(sub)   

pos_num                     = size(sub.xyz,1);

sub_origin                  = utils.get_GridOrigin(sub);
origined_xyz                = (sub.xyz-repmat(sub_origin,pos_num,1));

sub_basis                   = diag(max(sub.xyz)-min(sub.xyz));
grid_basis                  = diag([25,50,0]);
align_coeffs                = sub_basis\grid_basis;

aligned_xyz                 = align_coeffs*origined_xyz';


sub.aligned_xyz             = aligned_xyz';
%M1_position_xyz             = aligned_xyz'+repmat(utils.get_M1_xyz(),pos_num,1);
%sub.aligned_xyz             = M1_position_xyz;

end