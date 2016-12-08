function MEP = calculate_MEP(rawAmp,rawLAT)

whichalgorithm = 'Rossini'; %alternative: 'MillsNithi'

% Threshold MEPs
if strcmpi(whichalgorithm,'Rossini')
% Rossini-Rothwell: 
amp_threshold               = 50;
MEP                         =  single(rawAmp>amp_threshold);  
elseif strcmpi(whichalgorithm,'MillsNithi')
% Mills-Nithi
amp_flag                    = rawAmp>20;
lat_lowflag                 = rawLAT>17;
lat_hiflag                  = rawLAT<30;
MEP                         = single(amp_flag & lat_lowflag & lat_hiflag);  
else error('No correct algorithm specified');





end