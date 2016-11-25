function output = perform_weighting(sub,quad_weigths,threshed_weights)
    
    output = struct();
    
    output.quad_AMP      = quad_weigths*sub.amp;
    output.quad_MEP      = quad_weigths*sub.mep;     
    output.quad_LAT      = quad_weigths*sub.lat;    
    
    output.th_AMP      = threshed_weights*sub.amp;
    output.th_MEP      = threshed_weights*sub.mep;     
    output.th_LAT      = threshed_weights*sub.lat;    

end