function xyz = get_DesignGrid()
% Define Grid based on planned grid-points
% Step distance on rectangular grid therefore 3.5355 
% pdist2([0,0],[3.5355,3.5355]) = 5; sqrt(25/2) = 3.5355 -> 
STEPSIZE            = single(sqrt(25/2));
xyz                 = utils.get_M1_xyz;
do_center           = @(x)x-mean(x);
[X,Y,Z]             = meshgrid(do_center(0:STEPSIZE:6*STEPSIZE),do_center(0:STEPSIZE:14*STEPSIZE),0);
X                   = X+xyz(1);
Y                   = Y+xyz(2);
Z                   = Z+xyz(3);
do_structure        = @(x)reshape(x,1,[])';
xyz                 = cat(2,do_structure(X),do_structure(Y),do_structure(Z));

end