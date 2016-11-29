function plot_gridmodel(V,CLIM)

if nargin<2,
    CLIM                            = 2; %abs of upper and lower limit
end

% Prepare Data
V           = utils.vec2mesh(V);
Grid        = utils.get_DesignGrid;
X           = utils.vec2mesh(Grid(:,1));
Y           = utils.vec2mesh(Grid(:,2));
Z           = utils.vec2mesh(Grid(:,3));

% Bind Value to upper/lower bounds 

absTriValue                     = abs(V);
sgnTriValue                     = sign(V);    
absTriValue(absTriValue>CLIM)   = CLIM;    
BV                              = (absTriValue.*sgnTriValue)./CLIM;

% Interpolate
Xq          = linspace(min(Grid(:,1)),max(Grid(:,1)),150);
Yq          = linspace(min(Grid(:,2)),max(Grid(:,2)),70);
[Xq,Yq]     = meshgrid(Xq,Yq);
%BV          = interp2(X,Y,BV,Xq,Yq,'nearest');
%V           = interp2(X,Y,V,Xq,Yq,'nearest');
BV          = interp2(X,Y,BV,Xq,Yq,'cubic');
V           = interp2(X,Y,V,Xq,Yq,'cubic');
X           = Xq;
Y           = Yq;

% Plotting data
figure
hold on
h_c = contour(X,Y,BV,linspace(-CLIM,CLIM,40),'fill','on');
if any(any(abs(V)>1.30))   
    PV = double(abs(V)<1.301);
    h_s = contour(X,Y,PV,[0 .99],'fill','off','linecolor',[.5 .5 .5],'linewidth',2);        
end



caxis([-CLIM CLIM])
xlabel('X-Axis in mm')
ylabel('Y-Axis in mm')

% Configurate Colorbar
hnd_cb                      = colorbar;
[ytick,yticklab]            = make_yticks(CLIM);
set(hnd_cb,'YLIM',[-CLIM CLIM],'YTICK',ytick,'YTICKLABEL',yticklab)

% Plot Position of HotSpot according to Design (1 cm posterior to Grid Origin
hold on
xyz     = utils.get_DesignGridOrigin+[0,-10,0];
plot3(xyz(:,1),xyz(:,2),xyz(:,3),'ko','markerfacecolor','r')


end

function [ytick,yticklab] = make_yticks(CLIM)

    ytick                       = [-CLIM:-2,-1.301,0,1.301,2:CLIM];
    yticklab                    = {};
    for i_yt = 1:length(ytick)
        yticklab{i_yt} = sprintf('%.3g ',10^-abs(ytick(i_yt)));
    end
end


