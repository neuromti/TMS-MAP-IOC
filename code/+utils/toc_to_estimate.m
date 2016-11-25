function toc_to_estimate(tm,num_rep)

breakDown = @(x,scaler)[floor(x*scaler),(x*scaler)-floor(x*scaler)];

% Converted to weeks      
divisor     = (60*60*24*7);
fctr        = num_rep/divisor;
time_weeks  = (tm*fctr);

show_time       = [];
show_time(1,:)  = breakDown(time_weeks,1); %weeks
unit_scaler     = [1,7,24,60,60];
for unit_idx = 2:5, %days to seconds
show_time(unit_idx,:)  = breakDown(show_time(unit_idx-1,2),unit_scaler(unit_idx));
end

% Display estimated result   
fprintf('Will probably run for %0.0f weeks %0.0f days %0.0f hours %0.0f minutes and %0.0f seconds! \n',show_time(:,1)')

end