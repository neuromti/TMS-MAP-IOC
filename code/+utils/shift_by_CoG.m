function    sub = shift_by_CoG(sub)   

pos_num                     = size(sub.xyz,1);
%literature_based_M1         = [-37, -31, 71];
literature_based_M1         = [-38, -19 70]; %based on own data(average of all subjects CoG);
repeated_M1                 = repmat(literature_based_M1,pos_num,1);
originM1_xyz                = sub.xyz-repmat(sub.CoG,pos_num,1);
sub.shifted_xyz             = originM1_xyz+repeated_M1;

coeffs                      = pca(originM1_xyz);
%coeffs                      = sortrows(coeffs',-1)';
origin_realigned_xyz        = (coeffs*originM1_xyz')';
sub.aligned_xyz             = origin_realigned_xyz+repeated_M1;