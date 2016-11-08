%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [output]= filt_freqz(input,lower,upper,Fs,pass,extension,order)
% if not specified, Fs=1000;  
% if extension ~=0, symetric padarray
% checks automatically for bandpass, low or high
% written by M.Vukelic and R.Bauer for CIN AG NPT, Tuebingen 
% Version: 08.06.2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output]= filt_freqz(input,lower,upper,Fs,pass,extension,order)

if nargin < 4, Fs = 1000; end
if nargin < 5, pass = 2; end
if nargin <6, extension = 0; end
if nargin < 7, order = 0; end
if size(input,1) > size(input,2), input = input'; end % turn row to column

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% define parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

input   = double(input);
[C L]   = size(input);
output  = zeros(C,L);

if ~lower, lower = []; end
if ~upper, upper = []; end

fir_band = [2*(lower/Fs),2*(upper/Fs)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


if lower,
    if upper,    
        if ~order, order = max([3*ceil(Fs./[upper,lower]), 15]);  end
        if upper>lower,
            order = ceil(order/10)*10;
            b = fir1 (order,fir_band,'bandpass');
        elseif lower>upper,
            order = ceil(order/10)*10;
            b = fir1 (order,fliplr(fir_band),'stop');    
        else
            error('Strange Frequencies');
        end
    else
        if ~order, order = max([3*ceil(Fs./lower), 15]);  end
        order = ceil(order/10)*10;
        b = fir1 (order,fir_band,'high');   
    end
elseif upper,
    if ~order, order = max([3*ceil(Fs./upper), 15]);  end
    order = ceil(order/10)*10;
    b = fir1 (order,fir_band,'low');       
else
    error('No Frequencies');
end

if extension==0, 
    input = padarray(input,[0 ceil(order*1.5)],0,'both'); 
elseif extension == 2,
    input = padarray(input,[0 ceil(order*1.5)],'symmetric','both'); 
end
temp2 = zeros(size(input));
if pass == 2,
    for chan=1:C, temp2(chan,:)=filtfilt(b,1,input(chan,:)); end
elseif pass == 1,
     for chan=1:C, temp2(chan,:)=filter(b,1,input(chan,:)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% generate output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

output(1:C,1:L) = temp2(1:C,ceil(order*1.5)+1:ceil(order*1.5)+L);

end