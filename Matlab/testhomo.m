clc;
clear;
M=load('Model.txt');
m1=load('data1.txt');
m2=load('data2.txt');
m3=load('data3.txt');
M1=[M(:,1:2) ; M(:,3:4) ; M(:,5:6) ; M(:,7:8)];
m1=[m1(:,1:2) ; m1(:,3:4) ; m1(:,5:6) ; m1(:,7:8)];
m2=[m2(:,1:2) ; m2(:,3:4) ; m2(:,5:6) ; m2(:,7:8)];
m3=[m3(:,1:2) ; m3(:,3:4) ; m3(:,5:6) ; m3(:,7:8)];
mm(:,:,1)=m1';
mm(:,:,2)=m2';
mm(:,:,3)=m3';

[rows,npts]=size(M1');%npts 为列数
    matrixone=ones(1,npts);% 1矩阵
    M2=[M1';matrixone];%增加一行 1 1
    num=size(mm,3)
    for i=1:num
        mm(3,:,i)=matrixone; 
    end
x1=M2;
x2=mm;

Npts = length(x1);
    A = zeros(3*Npts,9);
    
    O = [0 0 0];
    for n = 1:Npts
	X = x1(:,n)';%定义 
	x = x2(1,n); y = x2(2,n); w = x2(3,n);
	A(3*n-2,:) = [  O  -w*X  y*X];
	A(3*n-1,:) = [ w*X   O  -x*X];
	A(3*n  ,:) = [-y*X  x*X   O ]
    end
    
    [U,D,V] = svd(A)
    % Ax=b  x=A\b;
    % Extract homography单应性矩阵
    H = reshape(V(:,9),3,3)';
    H=H/H(3,3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Maximun likelihood estimation for the H最大似然估计
    % using the function(10), P7
    %options = optimset('LargeScale','off','LevenbergMarquardt','on');
    %[x,resnorm,residual,exitflag,output]  = lsqnonlin( @simon_H, reshape(H,1,9) , [],[],options,mm, M2);
   %H=reshape(x,3,3);
    H1=H/H(3,3);
     V=[];
        v12(:,:)=[H(1,1)*H(2,1), H(1,1)*H(2,2)+H(1,2)*H(2,1), H(1,2)*H(2,2), H(1,3)*H(2,1)+H(1,1)*H(2,3), H(1,3)*H(2,2)+H(1,2)*H(2,3), H(1,3)*H(2,3)];
        v11(:,:)=[H(1,1)*H(1,1), H(1,1)*H(1,2)+H(1,2)*H(1,1), H(1,2)*H(1,2), H(1,3)*H(1,1)+H(1,1)*H(1,3), H(1,3)*H(1,2)+H(1,2)*H(1,3), H(1,3)*H(1,3)];
        v22(:,:)=[H(2,1)*H(2,1), H(2,1)*H(2,2)+H(2,2)*H(2,1), H(2,2)*H(2,2), H(2,3)*H(2,1)+H(2,1)*H(2,3), H(2,3)*H(2,2)+H(2,2)*H(2,3), H(2,3)*H(2,3)];
        V=[V;v12(:,:);v11(:,:)-v22(:,:)]
    
    k=V'*V;       
    [u,v,d]=svd(k);%奇异值分解[u,s,v]=svd(A),使得A=USV'
    [e,d2]=eig(k);%Eigenvector特征向量 [V,D]=eig(A)使得 AV=VD，D是特征值对角阵,V是特征向量阵
    b=d(:,6);%b就是论文作中B
    v0=(b(2)*b(4)-b(1)*b(5))/(b(1)*b(3)-b(2)^2);
    s=b(6)-(b(4)^2+v0*(b(2)*b(4)-b(1)*b(5)))/b(1);
    alpha_u=sqrt(s/b(1));
    alpha_v=sqrt(s*b(1)/(b(1)*b(3)-b(2)^2));
    skewness=-b(2)*alpha_u*alpha_u*alpha_v/s;
    u0=skewness*v0/alpha_u-b(4)*alpha_u*alpha_u/s;
    A=[alpha_u skewness u0
        0      alpha_v  v0
        0      0        1]
    B=rotate(pi/4,pi/6,pi/3);
    det(B);inv(B)*B
           