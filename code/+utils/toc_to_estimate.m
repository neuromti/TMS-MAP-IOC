function toc_to_estimate(tm,num_rep)
    % full time estimate in s
    ntm_in_s = tm*num_rep;
    rem_in_s = rem(ntm_in_s,60);
    % converted to min
    ntm_in_m = floor(ntm_in_s/60);
    rem_in_m = rem(ntm_in_m,60);
    % converted to h
    ntm_in_h = floor(ntm_in_m/60);
    rem_in_h = rem(ntm_in_h,60);
    % converted to d
    ntm_in_d = floor(ntm_in_h/24);
    rem_in_d = rem(ntm_in_d,24);
    
    fprintf('Will probably run for %i days %i hours %i minutes and %0.0f seconds! \n',rem_in_d,rem_in_h,rem_in_m,rem_in_s)
end