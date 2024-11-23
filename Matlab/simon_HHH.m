%
% ***********************************************************************************
% *******          A Flexible New Technique for Camera Calibration            *******
% ***********************************************************************************
%                            7/2004    
%                            
%
% REF:	   "A Flexible New Technique for Camera Calibration"
%           - Zhengyou Zhang 
%           - Microsoft Research 
%
function f = simon_HHH(params, m, M)
% unpack the params
  num=size(m,3);
  R=[];
  for i=1:num
      R_new=params( [(i-1)*6+1 : (i-1)*6+6] );
      Q1=R_new(1);
      Q2=R_new(2);
      Q3=R_new(3);
      TL=R_new([4:6])';
      RL=[cos(Q2)*cos(Q1)   sin(Q2)*cos(Q1)   -sin(Q1) ; -sin(Q2)*cos(Q3)+cos(Q2)*sin(Q1)*sin(Q3)    cos(Q2)*cos(Q3)+sin(Q2)*sin(Q1)*sin(Q3)  cos(Q1)*sin(Q3) ; sin(Q2)*sin(Q3)+cos(Q2)*sin(Q1)*cos(Q3)    -cos(Q2)*sin(Q3)+sin(Q2)*sin(Q1)*cos(Q3)  cos(Q1)*cos(Q3)];
      RT=[RL(:,1:2) , TL];
      R=[R;RT];
  end
  k1=params(num*6+1);
  k2=params(num*6+2);
  A=[params(num*6+3) params(num*6+4) params(num*6+5); 0 params(num*6+6) params(num*6+7); 0,0,1];
  u0=A(1,3);
  v0=A(2,3);
D=[];
d=[];
npts=size(m,2);
for flag=1:num
    RT=R([(flag-1)*3+1 : (flag-1)*3+3],:);
    XY=RT*M;
    UV=A*XY;
    UV=[UV(1,:)./UV(3,:); UV(2,:)./UV(3,:); UV(3,:)./UV(3,:)];
    XY=[XY(1,:)./XY(3,:); XY(2,:)./XY(3,:); XY(3,:)./XY(3,:)];
    D=[D; UV(1,:)+((UV(1,:)-u0).*( (XY(1,:)).^2 + (XY(2,:)).^2 ))*k1 + ((UV(1,:)-u0).*( (XY(1,:)).^2 + (XY(2,:)).^2 ).^2)*k2 ; UV(2,:) + ((UV(2,:)-v0).*( (XY(1,:)).^2 + (XY(2,:)).^2 ))*k1 + ((UV(2,:)-v0).*( (XY(1,:)).^2 + (XY(2,:)).^2 ).^2)*k2 ];
    d=[d; m(1,:,flag); m(2,:,flag)];
end
f=d-D;











