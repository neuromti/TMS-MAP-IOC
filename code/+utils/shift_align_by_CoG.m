function  sub = shift_align_by_CoG(sub)   

pos_num                     = size(sub.xyz,1);
M1_xyz                      = utils.get_M1_xyz();
repeated_M1                 = repmat(M1_xyz,pos_num,1);
originM1_xyz                = sub.xyz-repmat(sub.CoG,pos_num,1);

coeffs                      = pca(originM1_xyz);
%coeffs                     = sortrows(coeffs',-1)';
origin_realigned_xyz        = (coeffs*originM1_xyz')';
sub.shift_aligned_xyz       = origin_realigned_xyz+repeated_M1;