function R=Rotate(a,b,y)
Rx=[1   0     0
    0 cos(a) -sin(a)
    0 sin(a) cos(a)]
Ry=[ cos(b) 0 sin(b)
    0    1   0
    -sin(b) 0 cos(b)]
Rz=[cos(y) -sin(y) 0
    sin(y) cos(y) 0
    0   0    1]
R=Rx*Ry*Rz

