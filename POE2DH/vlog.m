function xi = vlog( g )
%vlog Logarithm mapping from a homogeneous transformation to a twist
%   
%   xi = vlog (g)
%   g:      Homogeneous transformation, 4 x 4
%   xi:     twist, 6 x 1

DELTA=10^(-12);

fai=rotationTheta(g);

w=vlogR(g(1:3,1:3));

if abs(fai)<DELTA
    
    p=g(1:3,4);
    
else
    
    p=(eye(3)-0.5*skew(w)+(2*sin(fai)-fai*(1+cos(fai)))/(2*fai*fai*sin(fai))*skew(w)*skew(w))*g(1:3,4);

end

xi=[w;p];

end

