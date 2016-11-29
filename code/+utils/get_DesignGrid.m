function xyz = get_DesignGrid()

% Define Grid based on planned grid-points
% Step distance on rectangular grid therefore 3.5355 
% pdist2([0,0],[3.5355,3.5355]) = 5; sqrt(25/2) = 3.5355 -> 
STEPNUM             = [7,15];
STEPSIZE            = single(sqrt(25/2));  

%xyz                 = utils.get_M1();
xyz                 = utils. get_DesignGridOrigin();

do_center           = @(x)x-mean(x);

StepsX              = do_center(linspace(0,(STEPNUM(1)-1)*STEPSIZE,STEPNUM(1)));
StepsY              = do_center(linspace(0,(STEPNUM(2)-1)*STEPSIZE,STEPNUM(2)));

[X,Y,Z]             = meshgrid(StepsX,StepsY,0);
X                   = X+xyz(1);
Y                   = Y+xyz(2);
Z                   = Z+xyz(3);
do_structure        = @(x)reshape(x,1,[])';
xyz                 = cat(2,do_structure(X),do_structure(Y),do_structure(Z));

end