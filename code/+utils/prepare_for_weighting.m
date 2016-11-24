function sub = prepare_for_weighting(mapping)    
amp_threshold               = 50;

remove_flag                 = mapping.amp<0 | mapping.lat <0;

sub.xyz                     = mapping.xyz(~remove_flag,:);
sub.amp                     = mapping.amp(~remove_flag);    
sub.lat                     = mapping.lat(~remove_flag);
sub.mep                     = double(sub.amp>amp_threshold);        
sub.CoG                     = (sub.amp./sum(sub.amp))'*sub.xyz;
%center the latency 
ispositive_flag             = sub.lat>0;
sub.lat(sub.lat>0)          = sub.lat(ispositive_flag)-mean(sub.lat(ispositive_flag));

    
    