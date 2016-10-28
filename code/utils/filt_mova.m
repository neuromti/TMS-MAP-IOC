%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function output = filt_mova (input,window,mode);
%
% if nargin < 2, window = 100; end
% if nargin < 3, mode = 1; end  
% if mode == 1, nanmean
% if mode == 2, nanmedian
% if mode == 3, root nanmean square
% if mode == 4, nanmean slope
% if mode == 5, nanmean absolute slope
% if mode == 6, nanmean baselined
% if mode == 7, erosion
% if mode == 8, max filter
% if mode == 9, gaussian filter
% if mode == 10, exponential filter
% if mode == 11, correlation
% written by R.Bauer for CIN AG NPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = filt_mova (input,window,mode,pad)

if nargin < 2, window = length(input)./10; end  
if nargin < 3, mode = 1; end  
if nargin < 4, pad = 0; end

if window == 0 || window == 1, output = input; return; end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pad~=0,
    O = size(input,2);
    if pad == 1,
        input   = padarray(input,[0,window],'replicate'); 
    elseif pad == 2,
        input   = padarray(input,[0,window],'symmetric'); 
    elseif pad == 3,
        input   = padarray(input,[0,window],'circular'); 
    end        
end

output = zeros(size(input,1), size(input,2)); %predefine for speed

L = size(input,2);

for n= 1:size(input,2),
    
a = n-(floor(window/2));    
b = n+(ceil(window/2));     

if a < 1, a=1; end 
if b > L, b=L; end;

for c= 1:size(input,1),
    if mode == 1, output(c,n) = nanmean(input(c,a:b)); end
    if mode == 2, output(c,n) = nanmedian(input(c,a:b)); end
    if mode == 3, output(c,n) = (nanmean(input(c,a:b).^2))^0.5; end
    if mode == 4, output(c,n) = nanmean(gradient(input(c,a:b))); end
    if mode == 5, output(c,n) = nanmean(abs(gradient(input(c,a:b)))); end
    if mode == 6, output(c,n) = nanmean(baseline(input(c,a:b),1)); end
    if mode == 7, output(c,n) = min(input(c,a:b)); end
    if mode == 8, output(c,n) = max(input(c,a:b)); end
    if mode == 9, output(c,n) = 2.*nanmean(gausswin(length(a:b))'.*(input(c,a:b))); end    
end

if mode == 11, output = zeros(size(input,1), size(input,1), size(input,2)); output(:,:,n) = corr(input(:,a:b)'); end
if mode == 10, output = conv(input,(1./exp(1:window))./sum(1./exp(1:window)),'same'); end
end

if pad~=0,
    output = output(:,window+1:window+O); 
end
    
end