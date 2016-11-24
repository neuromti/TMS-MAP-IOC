function    sub = shift_by_CoG(sub)
       
literature_based_M1         = [-37, -31, 71];
pos_num                     = size(sub.xyz,1);
sub.shifted_xyz             = sub.xyz-repmat(sub.CoG,pos_num,1)+repmat(literature_based_M1,pos_num,1);
