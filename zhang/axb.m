clc;
clear;
x1=[0   0.5  0.5  
   -0.5 -0.5  0   
   1    1     1  ];
x2=[ 63.43   92.46   91.8    
     405.5  407.45  438.65   
     1       1       1          ];
Npts = length(x1);
    A = zeros(3*Npts,9);
    
    O = [0 0 0];
    for n = 1:Npts
	X = x1(:,n)';%定义 
	x = x2(1,n); y = x2(2,n); w = x2(3,n);
	A(3*n-2,:) = [  O  -w*X  y*X];
	A(3*n-1,:) = [ w*X   O  -x*X];
	A(3*n  ,:) = [-y*X  x*X   O ]
    A
    end
    [U,D,V] = svd(A);
    U
    D
    V
    
    % Ax=b  x=A\b;
    % Extract homography单应性矩阵
    H1 = reshape(V(:,9),3,3)';
    H=H1/H1(3,3)
    