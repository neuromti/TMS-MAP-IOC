function cenLat = center_latency(latency)

cenLat                      = latency;
ispositive_flag             = latency>0;
cenLat(latency>0)           = latency(ispositive_flag)-mean(latency(ispositive_flag));

end