%
% ***********************************************************************************
% *******          A Flexible New Technique for Camera Calibration            *******
% ***********************************************************************************
%                            7/2004    Simon Wan
%                            //2006-03-04 如有疑问：simonwan1980@gmail.com (因为已从哈工大毕业，此地址已作废simonwan1980@hit.edu.cn)
%
% REF:	   "A Flexible New Technique for Camera Calibration"
%           - Zhengyou Zhang 
%           - Microsoft Research 
%
% HOMOGRAPHY2D - computes 2D homography
%
% Usage:   H = homography2d(x1, x2)
%          H = homography2d(x)
%
% Arguments:
%          x1  - 3xN set of homogeneous points
%          x2  - 3xN set of homogeneous points such that x1<->x2
%         
%           x  - If a single argument is supplied it is assumed that it
%                is in the form x = [x1; x2]
% Returns:
%          H - the 3x3 homography such that x2 = H*x1
%
% This code follows the normalised direct linear transformation 
% algorithm given by Hartley and Zisserman "Multiple View Geometry in
% Computer Vision" p92.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 本子程序来源于下：
% Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/~pk
%
% May 2003  - Original version.
% Feb 2004  - single argument allowed for to enable use with ransac.
%H=A*[r1,r2,t];
function H = homography2d(varargin)
    
    [x1, x2] = checkargs(varargin(:));% varargin"变长度输入宗量列表"varargin本身是个元胞数组
    M=x1;                             % varargout"变长度输出宗量列表"
    m=x2;
    % Attempt to normalise（ 规格化）each set of points so that the origin 
    % is at centroid （质心）and mean distance from origin is sqrt(2).（因为是正方形）
    [x1, T1] = normalise2dpts(x1);
    [x2, T2] = normalise2dpts(x2);
    
    % Note that it may have not been possible to normalise
    % the points if one was at infinity so the following does not
    % assume that scale parameter w = 1.
    % Estimation of the H between the model plane and its image, P18建立单应性矩阵
    Npts = length(x1);
    A = zeros(3*Npts,9);%A为超定方程
    
    O = [0 0 0];
    for n = 1:Npts
	X = x1(:,n)';%定义 
	x = x2(1,n);y = x2(2,n); w = x2(3,n);
	A(3*n-2,:) = [  O  -w*X  y*X];
	A(3*n-1,:) = [ w*X   O  -x*X];
	A(3*n  ,:) = [-y*X  x*X   O ];
    end
    
    [U,D,V] = svd(A);
    % Ax=b  x=A\b;
    % Extract homography单应性矩阵
    H1 = reshape(V(:,9),3,3)'
           
    % Denormalize反向规格化,
    H2= T2\H1*T1;
    H=H2/H2(3,3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Maximun likelihood estimation for the H最大似然估计
    % using the function(10), P7
    options = optimset('LargeScale','off','LevenbergMarquardt','on');
    [x,resnorm,residual,exitflag,output]  = lsqnonlin( @simon_H, reshape(H,1,9) , [],[],options,m, M);
    H=reshape(x,3,3);
    H=H/H(3,3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% Function to check argument values and set defaults

function [x1, x2] = checkargs(arg);
    
    if length(arg) == 2
	x1 = arg{1};
	x2 = arg{2};
	if ~all(size(x1)==size(x2))
	    error('x1 and x2 must have the same size');
	elseif size(x1,1) ~= 3
	    error('x1 and x2 must be 3xN');
	end
	
    elseif length(arg) == 1
	if size(arg{1},1) ~= 6
	    error('Single argument x must be 6xN');
	else
	    x1 = arg{1}(1:3,:);
	    x2 = arg{1}(4:6,:);
	end
    else
	error('Wrong number of arguments supplied');
    end
    