function  sub = normalize_by_HsAnt(sub,setup)   

pos_num                     = size(sub.xyz,1);

if any(isnan(sub.ANT)) || any(isnan(sub.HS))
    sub.aligned_xyz             = NaN(pos_num,3);
    warning('Jumped a Subject because no proper Anterior or HotSpot was recorded!');
    return;
end

group_M1                    = utils.get_M1_xyz();
sub_M1                      = utils.calculate_CoG(sub);

sub_origin                  = nanmean(cat(1,sub.ANT,sub.HS));
sub_basis                   = sub_origin-sub.ANT;
group_origin                = nanmean(nanmean(cat(3,setup.ANT,setup.HS),3),1);
group_basis                 = group_origin-nanmean(setup.ANT);

align_coeffs                = diag(sub_basis)\diag(group_basis); % find the rescaling factor to align (based on subejcts and groups basis)

origined_xyz                = sub.xyz-repmat(sub_origin,pos_num,1); %translate to origin for all of subjects grid points
aligned_xyz                 = align_coeffs*origined_xyz'; % rescale with factor to align
reposition_xyz              = aligned_xyz'+repmat(group_origin,pos_num,1);
M1_position_xyz             = reposition_xyz+repmat((group_M1-sub_M1),pos_num,1);

sub.aligned_xyz             = zscore(M1_position_xyz);