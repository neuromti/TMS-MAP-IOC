function    sub = shift_by_CoG(sub)   

pos_num                     = size(sub.xyz,1);
M1_xyz                      = utils.get_M1_xyz();
repeated_M1                 = repmat(M1_xyz,pos_num,1);
originM1_xyz                = sub.xyz-repmat(sub.CoG,pos_num,1);
sub.shifted_xyz             = originM1_xyz+repeated_M1;
