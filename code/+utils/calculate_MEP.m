function MEP = calculate_MEP(rawAmp,rawLatency)

if nargin<2, rawLatency = false; end

amp_threshold               = 50;
MEP                         =  single(rawAmp>amp_threshold);  

end