function plot_landmark2grid(name,coordsys)    
    
if nargin <2, coordsys = 'MNI'; end
if nargin <1, name = {'S1','M1','PMv','PMd','SMAproper','preSMA'}; end

for lm_i = 1 : length(name)  
    [xyz,name] = utils.get_landmark(name,coordsys);
    plot3(xyz(lm_i,1),xyz(lm_i,2),xyz(lm_i,3),'ko','markerfacecolor','r')
end

end