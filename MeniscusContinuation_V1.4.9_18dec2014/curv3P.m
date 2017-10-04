function Ka = curv3P(x1,y1,x2,y2,x3,y3)
% inverse of the radius of a circle going through 3 points
Ka = 2.*(x1*y2-(x1)*y3-(x2)*y1+x2*y3+x3*y1-(x3)*y2)...
    /sqrt(((x2-x1)^2+(y2-y1)^2)*((x2-x3)^2+(y2-y3)^2)*((x3-x1)^2+(y3-y1)^2) );

