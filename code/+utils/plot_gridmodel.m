function plot_gridmodel(V,CLIM)

if nargin<2,
    CLIM                            = 2; %abs of upper and lower limit
end

% Prepare Data
V           = utils.vec2mesh(V);
Grid        = utils.get_DesignGrid;
X           = utils.vec2mesh(Grid(:,1));
Y           = utils.vec2mesh(Grid(:,2));
% Z           = utils.vec2mesh(Grid(:,3));

% Bind Value to upper/lower bounds 

absTriValue                     = abs(V);
sgnTriValue                     = sign(V);    
absTriValue(absTriValue>CLIM)   = CLIM;    
BV                              = (absTriValue.*sgnTriValue)./CLIM;

% Interpolate
Xq          = linspace(min(Grid(:,1)),max(Grid(:,1)),14);
Yq          = min(Grid(:,2)):mean(diff(Xq)):max(Grid(:,2));

[Xq,Yq]     = meshgrid(Xq,Yq);
iBV         = interp2(X,Y,BV,Xq,Yq,'linear');
iV          = interp2(X,Y,V,Xq,Yq,'linear');
X           = Xq;
Y           = Yq;

% Plotting data
figure
hold on
contour(X,Y,iBV,linspace(-CLIM,CLIM,40),'fill','on');
if any(any(abs(iV)>1.30))       
    contour(X,Y,abs(iV),[0 1.301],'fill','off','linecolor',[.5 .5 .5],'linewidth',2);        
end

caxis([-CLIM CLIM])
xlabel('X-Axis in mm')
ylabel('Y-Axis in mm')

% Configurate Colorbar
hnd_cb                      = colorbar;
[ytick,yticklab]            = make_yticks(CLIM);
set(hnd_cb,'YLIM',[-CLIM CLIM],'YTICK',ytick,'YTICKLABEL',yticklab)

% Plot Position of HotSpot according to Design (1 cm posterior to Grid Origin
target_flat = false;
if target_flat
    hold on
    HS          = utils.get_DesignGridOrigin+[0,-10,0];
    plot3(HS(:,1),HS(:,2),HS(:,3),'ko','markerfacecolor','r')

    ANT         = utils.get_GroupAnt();
    plot3(ANT(:,1),ANT(:,2),ANT(:,3),'ko','markerfacecolor','b')
end

end

function [ytick,yticklab] = make_yticks(CLIM)
    ytick                       = [-CLIM:-2,-1.301,0,1.301,2:CLIM];
    yticklab                    = {};
    for i_yt = 1:length(ytick)
        yticklab{i_yt} = sprintf('%.3g ',10^-abs(ytick(i_yt)));
    end
end