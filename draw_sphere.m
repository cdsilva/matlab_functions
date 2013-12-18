function draw_sphere(x,y,z,color,r)
%%function draws sphere at position (x,y,z) with radius specified by r

[X,Y,Z]=sphere;

surf(r*X+x,r*Y+y,r*Z+z,'facecolor',color)
