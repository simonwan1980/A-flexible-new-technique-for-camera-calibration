% function Zhang(M,m)
%
% ***********************************************************************************
% *******          A Flexible New Technique for Camera Calibration            *******
% ***********************************************************************************
%                            7/2004    Simon Wan
%                            //2006-03-04 如有疑问：simonwan1980@gmail.com (因为已从哈工大毕业，此地址已作废simonwan1980@hit.edu.cn)
%
% Note:    M:2*N  m:2*N
% M        point on the model plane, when using M=[X,Y]' ---> M=[X,Y,1]'
% m        M's image, when using                m=[u,v]' ---> m=[u,v,1]' , so that
%          s*m = H*M , with H=A*[r1,r2,t];                  (2)
% H        homography matrix
%
% REF:	   "A Flexible New Technique for Camera Calibration"
%           - Zhengyou Zhang 
%           - Microsoft Research 
%
function Zhang(M,m)

%  M=[X,Y]' ---> M=[X,Y,1]'  ;   m=[u,v]' ---> m=[u,v,1]' 
    [rows,npts]=size(M);
    matrixone=ones(1,npts);% 1矩阵
    M=[M;matrixone];%%3*256   %  M=[X,Y]' ---> M=[X,Y,1]' 
    num=size(m,3);%%m矩阵大小的第三个元素5
    for i=1:num
        m(3,:,i)=matrixone; %m=[u,v]' ---> m=[u,v,1]' 
    end
% Estimate the H
 %H=A*[r1,r2,t];
    for i=1:num
        H(:,:,i)=homography2d1(M,m(:,:,i))';%%%调用函数homography2d1.m
    end
% solve the intrinsic parameters matrix A
% A=[alpha_u skewness u0
%    0       alpha_v  v0
%    0       0        1]
% see Appendix B "Extraction of the Intrisic Parameters from Matrix B", P18
    V=[];
    for flag=1:num
        v12(:,:,flag)=[H(1,1,flag)*H(2,1,flag), H(1,1,flag)*H(2,2,flag)+H(1,2,flag)*H(2,1,flag), H(1,2,flag)*H(2,2,flag), H(1,3,flag)*H(2,1,flag)+H(1,1,flag)*H(2,3,flag), H(1,3,flag)*H(2,2,flag)+H(1,2,flag)*H(2,3,flag), H(1,3,flag)*H(2,3,flag)];
        v11(:,:,flag)=[H(1,1,flag)*H(1,1,flag), H(1,1,flag)*H(1,2,flag)+H(1,2,flag)*H(1,1,flag), H(1,2,flag)*H(1,2,flag), H(1,3,flag)*H(1,1,flag)+H(1,1,flag)*H(1,3,flag), H(1,3,flag)*H(1,2,flag)+H(1,2,flag)*H(1,3,flag), H(1,3,flag)*H(1,3,flag)];
        v22(:,:,flag)=[H(2,1,flag)*H(2,1,flag), H(2,1,flag)*H(2,2,flag)+H(2,2,flag)*H(2,1,flag), H(2,2,flag)*H(2,2,flag), H(2,3,flag)*H(2,1,flag)+H(2,1,flag)*H(2,3,flag), H(2,3,flag)*H(2,2,flag)+H(2,2,flag)*H(2,3,flag), H(2,3,flag)*H(2,3,flag)];
        V=[V;v12(:,:,flag);v11(:,:,flag)-v22(:,:,flag)];
    end
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
        0      0        1];
% solve k1 k1 and all the extrisic parameters, P6
    D=[];
    d=[];
    Rm=[];
    for flag=1:num
        s=(1/norm(inv(A)*H(1,:,flag)')+1/norm(inv(A)*H(2,:,flag)'))/2;
        rl1=s*inv(A)*H(1,:,flag)';%r1=s*inv(A)*h1
        rl2=s*inv(A)*H(2,:,flag)';%r2=s*inv(A)*h2
        rl3=cross(rl1,rl2);    %r3=r1Xr2;
                               %C = cross(A,B) returns the cross product of the vectors
                               %A and B.  That is, C = A x B.  A and B must be 3 element    vectors.
        RL=[rl1,rl2,rl3];
        %%%%%%%%%%%%%%%%%%%%
        % see Appendix C "Approximating a 3*3 matrix by a Rotation Matrix", P19
        [U,S,V] = svd(RL);
        RL=U*V';
        %%%%%%%%%%%%%%%%%%%%
        TL=s*inv(A)*H(3,:,flag)';%TL是位移矩阵t(外参)
        RT=[rl1,rl2,TL];%H=A[r1 r2 t]
        XY=RT*M;%M是model plane 点的坐标
        UV=A*XY;%sm=A[R t]M,UV是等式的右边
        UV=[UV(1,:)./UV(3,:); UV(2,:)./UV(3,:); UV(3,:)./UV(3,:)];
        XY=[XY(1,:)./XY(3,:); XY(2,:)./XY(3,:); XY(3,:)./XY(3,:)];
        for j=1:npts
            D=[D; ((UV(1,j)-u0)*( (XY(1,j))^2 + (XY(2,j))^2 )) , ((UV(1,j)-u0)*( (XY(1,j))^2 + (XY(2,j))^2 )^2) ; ((UV(2,j)-v0)*( (XY(1,j))^2 + (XY(2,j))^2 )) , ((UV(2,j)-v0)*( (XY(1,j))^2 + (XY(2,j))^2 )^2) ];
            d=[d; (m(1,j,flag)-UV(1,j)) ; (m(2,j,flag)-UV(2,j))];
        end
        r13=RL(1,3);
        r12=RL(1,2);
        r23=RL(2,3);
        Q1=-asin(r13);
        Q2=asin(r12/cos(Q1));%asin就是arcsin
        Q3=asin(r23/cos(Q1));
        [cos(Q2)*cos(Q1)   sin(Q2)*cos(Q1)   -sin(Q1) ; -sin(Q2)*cos(Q3)+cos(Q2)*sin(Q1)*sin(Q3)    cos(Q2)*cos(Q3)+sin(Q2)*sin(Q1)*sin(Q3)  cos(Q1)*sin(Q3) ; sin(Q2)*sin(Q3)+cos(Q2)*sin(Q1)*cos(Q3)    -cos(Q2)*sin(Q3)+sin(Q2)*sin(Q1)*cos(Q3)  cos(Q1)*cos(Q3)];
        R_new=[Q1,Q2,Q3,TL'];
        Rm=[Rm , R_new];
    end
% using function (13), P8
    k=inv(D'*D)*D'*d;
% Complete Maximun Likelihood Estimation, using function (14), P8
    para=[Rm,k(1),k(2),alpha_u,skewness,u0,alpha_v,v0];
    % optimset Create/alter OPTIM OPTIONS structure.
    % options = optimset('LargeScale','off','LevenbergMarquardt','on');
    options.Algorithm = 'levenberg-marquardt';
    %lsqnonlin :Solves non-linear least squares problems.最小二乘法问题
         %     Examples
         %    FUN can be specified using @:
         %    x = lsqnonlin(@myfun,[2 3 4])
         %  where MYFUN is a MATLAB function such as:
         %        function F = myfun(x)
         %  F = sin(x);
         %    FUN can also be an inline object:
         %        fun = inline('sin(3*x)')
         %  x = lsqnonlin(fun,[1 4]);
    [x,resnorm,residual,exitflag,output]  = lsqnonlin( @simon_HHH, para, [],[],options, m, M);
% display the result
    k1=x(num*6+1)
    k2=x(num*6+2)
    A=[x(num*6+3) x(num*6+4) x(num*6+5); 0 x(num*6+6) x(num*6+7); 0,0,1]




