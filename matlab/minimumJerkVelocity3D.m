function [Bx,By,Bz,B,Jx,Jy,Jz,J,Hx,Hy,Hz,H] = minimumJerkVelocity3D(t0,D,Ax,Ay,Az,t)
% MJxyz - evaluate a minimum jerk curve with seperate displacement for xyz
%
% [Bx,By,Bz,B,Jx,Jy,Jz,J,Hx,Hy,Hz,H] = MJxy(t0,D,Ax,Ay,Az,t)
%
% see Flash and Hogan (1985) for details on the minimum jerk equation
%
%    t0 = movement start time
%    D  = movement duration
%    Ax = displacement resulting from the movement (x)
%    Ay = displacement resulting from the movement (y)
%    Az = displacement resulting from the movement (y)
%
% The function is evaluated at times t
%
% The function also optionally returns the first-order and second-order
% partial derivatives, for use with optimization routines
%
% Bx, By, Bz and B are the x velocity, y velocity, z velocity and tangential velocities
% Jx, Jy, Jz and J are the gradients (partial derivatives) of the same quantities
% Hx, Hy, Hz and H are the Hessian (second-order partial derivatives)
%

% Jason Friedman, 2021
% www.curiousjason.com


% nt is normalized time (0 <= nt <= 1)
nt = (t(r)-t0)./D;
r = (nt>=0 & nt<=1);
Bx = zeros(1,numel(t));
By = zeros(1,numel(t));
Bz = zeros(1,numel(t));
B = zeros(1,numel(t));

Bx(r) = Ax/D * (-60 * nt(r).^3 + 30 * nt(r).^4 + 30 * nt(r).^2);
By(r) = Ay/D * (-60 * nt(r).^3 + 30 * nt(r).^4 + 30 * nt(r).^2);
Bz(r) = Az/D * (-60 * nt(r).^3 + 30 * nt(r).^4 + 30 * nt(r).^2);

A_tang = sqrt((Ax/D).^2 + (Ay/D).^2 + (Az/D).^2);
B(r) = A_tang * (-60 * nt(r).^3 + 30 * nt(r).^4 + 30 * nt(r).^2);

if nargout > 3
    Jx = zeros(5,numel(t));
    Jy = zeros(5,numel(t));
    Jz = zeros(5,numel(t));
    J = zeros(5,numel(t));
    
    Jx(1,r) = Ax.*(-(60.*(t(r)-t0))./D.^3+(180.*(t(r)-t0).^2)./D.^4-(120.*(t(r)-t0).^3)./D.^5);
    Jy(1,r) = Ay.*(-(60.*(t(r)-t0))./D.^3+(180.*(t(r)-t0).^2)./D.^4-(120.*(t(r)-t0).^3)./D.^5);
    Jz(1,r) = Az.*(-(60.*(t(r)-t0))./D.^3+(180.*(t(r)-t0).^2)./D.^4-(120.*(t(r)-t0).^3)./D.^5);

    Jx(2,r) = Ax.*(-(90.*(t(r)-t0).^2)./D.^4+(240.*(t(r)-t0).^3)./D.^5-(150.*(t(r)-t0).^4)./D.^6);
    Jy(2,r) = Ay.*(-(90.*(t(r)-t0).^2)./D.^4+(240.*(t(r)-t0).^3)./D.^5-(150.*(t(r)-t0).^4)./D.^6);
    Jz(2,r) = Az.*(-(90.*(t(r)-t0).^2)./D.^4+(240.*(t(r)-t0).^3)./D.^5-(150.*(t(r)-t0).^4)./D.^6);
    
    Jx(3,r) = (30.*(t(r)-t0).^2)./D.^3-(60.*(t(r)-t0).^3)./D.^4+(30.*(t(r)-t0).^4)./D.^5;
    
    Jy(4,r) = (30.*(t(r)-t0).^2)./D.^3-(60.*(t(r)-t0).^3)./D.^4+(30.*(t(r)-t0).^4)./D.^5;

    Jz(5,r) = (30.*(t(r)-t0).^2)./D.^3-(60.*(t(r)-t0).^3)./D.^4+(30.*(t(r)-t0).^4)./D.^5;

    
    J(1,r) = A_tang.*(-(60.*(t(r)-t0))./D.^2+(180.*(t(r)-t0).^2)./D.^3-(120.*(t(r)-t0).^3)./D.^4);
    J(2,r) = ((-(2.*Az.^2)./D.^3-(2.*Ay.^2)./D.^3-(2.*Ax.^2)./D.^3).*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(2.*A_tang)+(-(60.*(t(r)-t0).^2)./D.^3+(180.*(t(r)-t0).^3)./D.^4-(120.*(t(r)-t0).^4)./D.^5).*A_tang;
    J(3,r) = (Ax.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(A_tang.*D.^2);
    J(4,r) = (Ay.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(A_tang.*D.^2);
    J(5,r) = (Az.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(A_tang.*D.^2);

    if nargout >6
        Hx(1,1,:) = Ax.*(60./D.^3-(360.*(t(r)-t0))./D.^4+(360.*(t(r)-t0).^2)./D.^5);
        Hy(1,1,:) = Ay.*(60./D.^3-(360.*(t(r)-t0))./D.^4+(360.*(t(r)-t0).^2)./D.^5);
        Hz(1,1,:) = Az.*(60./D.^3-(360.*(t(r)-t0))./D.^4+(360.*(t(r)-t0).^2)./D.^5);
        H(1,1,:) = (60./D.^2-(360.*(t(r)-t0))./D.^3+(360.*(t(r)-t0).^2)./D.^4).*A_tang;
        
        Hx(2,1,:) = Ax.*((180.*(t(r)-t0))./D.^4-(720.*(t(r)-t0).^2)./D.^5+(600.*(t(r)-t0).^3)./D.^6);
        Hy(2,1,:) = Ay.*((180.*(t(r)-t0))./D.^4-(720.*(t(r)-t0).^2)./D.^5+(600.*(t(r)-t0).^3)./D.^6);
        Hz(2,1,:) = Az.*((180.*(t(r)-t0))./D.^4-(720.*(t(r)-t0).^2)./D.^5+(600.*(t(r)-t0).^3)./D.^6);
        H(2,1,:) = ((-(2.*Az.^2)./D.^3-(2.*Ay.^2)./D.^3-(2.*Ax.^2)./D.^3).*(-(60.*(t(r)-t0))./D.^2+(180.*(t(r)-t0).^2)./D.^3-(120.*(t(r)-t0).^3)./D.^4))./(2.*A_tang)+((120.*(t(r)-t0))./D.^3-(540.*(t(r)-t0).^2)./D.^4+(480.*(t(r)-t0).^3)./D.^5).*A_tang;
        
        Hx(3,1,:) = -(60.*(t(r)-t0))./D.^3+(180.*(t(r)-t0).^2)./D.^4-(120.*(t(r)-t0).^3)./D.^5;
        Hy(3,1,:) = zeros(size(t));
        Hz(3,1,:) = zeros(size(t));
        H(3,1,:) = (Ax.*(-(60.*(t(r)-t0))./D.^2+(180.*(t(r)-t0).^2)./D.^3-(120.*(t(r)-t0).^3)./D.^4))./(A_tang.*D.^2);
        
        Hx(4,1,:) = zeros(size(t));
        Hy(4,1,:) = -(60.*(t(r)-t0))./D.^3+(180.*(t(r)-t0).^2)./D.^4-(120.*(t(r)-t0).^3)./D.^5;
        Hz(4,1,:) = zeros(size(t));
        H(4,1,:) = (Ay.*(-(60.*(t(r)-t0))./D.^2+(180.*(t(r)-t0).^2)./D.^3-(120.*(t(r)-t0).^3)./D.^4))./(A_tang.*D.^2);
      
        Hx(5,1,:) = zeros(size(t));
        Hy(5,1,:) = zeros(size(t));
        Hz(5,1,:) = -(60.*(t(r)-t0))./D.^3+(180.*(t(r)-t0).^2)./D.^4-(120.*(t(r)-t0).^3)./D.^5;
        H(5,1,:) = (Az.*(-(60.*(t(r)-t0))./D.^2+(180.*(t(r)-t0).^2)./D.^3-(120.*(t(r)-t0).^3)./D.^4))./(A_tang.*D.^2);
       
        Hx(1,2,:) = Ax.*((180.*(t(r)-t0))./D.^4-(720.*(t(r)-t0).^2)./D.^5+(600.*(t(r)-t0).^3)./D.^6);    
        Hy(1,2,:) = Ay.*((180.*(t(r)-t0))./D.^4-(720.*(t(r)-t0).^2)./D.^5+(600.*(t(r)-t0).^3)./D.^6);
        Hz(1,2,:) = Az.*((180.*(t(r)-t0))./D.^4-(720.*(t(r)-t0).^2)./D.^5+(600.*(t(r)-t0).^3)./D.^6);
        H(1,2,:) = ((-(2.*Az.^2)./D.^3-(2.*Ay.^2)./D.^3-(2.*Ax.^2)./D.^3).*(-(60.*(t(r)-t0))./D.^2+(180.*(t(r)-t0).^2)./D.^3-(120.*(t(r)-t0).^3)./D.^4))./(2.*A_tang)+((120.*(t(r)-t0))./D.^3-(540.*(t(r)-t0).^2)./D.^4+(480.*(t(r)-t0).^3)./D.^5).*A_tang;
        
        Hx(2,2,:) = Ax.*((360.*(t(r)-t0).^2)./D.^5-(1200.*(t(r)-t0).^3)./D.^6+(900.*(t(r)-t0).^4)./D.^7);
        Hy(2,2,:) = Ay.*((360.*(t(r)-t0).^2)./D.^5-(1200.*(t(r)-t0).^3)./D.^6+(900.*(t(r)-t0).^4)./D.^7);
        Hz(2,2,:) = Az.*((360.*(t(r)-t0).^2)./D.^5-(1200.*(t(r)-t0).^3)./D.^6+(900.*(t(r)-t0).^4)./D.^7);
        H(2,2,:) = (((6.*Az.^2)./D.^4+(6.*Ay.^2)./D.^4+(6.*Ax.^2)./D.^4).*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(2.*A_tang)-((-(2.*Az.^2)./D.^3-(2.*Ay.^2)./D.^3-(2.*Ax.^2)./D.^3).^2.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(4.*(Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2))+((180.*(t(r)-t0).^2)./D.^4-(720.*(t(r)-t0).^3)./D.^5+(600.*(t(r)-t0).^4)./D.^6).*A_tang +((-(2.*Az.^2)./D.^3-(2.*Ay.^2)./D.^3-(2.*Ax.^2)./D.^3).*(-(60.*(t(r)-t0).^2)./D.^3+(180.*(t(r)-t0).^3)./D.^4-(120.*(t(r)-t0).^4)./D.^5))./A_tang;
        
        Hx(3,2,:) = -(90.*(t(r)-t0).^2)./D.^4+(240.*(t(r)-t0).^3)./D.^5-(150.*(t(r)-t0).^4)./D.^6;
        Hy(3,2,:) = zeros(size(t));
        Hz(3,2,:) = zeros(size(t));
        H(3,2,:) = -(Ax.*(-(2.*Az.^2)./D.^3-(2.*Ay.^2)./D.^3-(2.*Ax.^2)./D.^3).*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(2.*(Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2).*D.^2)+(Ax.*(-(60.*(t(r)-t0).^2)./D.^3+(180.*(t(r)-t0).^3)./D.^4-(120.*(t(r)-t0).^4)./D.^5))./(A_tang.*D.^2)-(2.*Ax.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(A_tang.*D.^3);
        
        Hx(4,2,:) = zeros(size(t));
        Hy(4,2,:) = -(90.*(t(r)-t0).^2)./D.^4+(240.*(t(r)-t0).^3)./D.^5-(150.*(t(r)-t0).^4)./D.^6;
        Hz(4,2,:) = zeros(size(t));
        H(4,2,:) = -(Ay.*(-(2.*Az.^2)./D.^3-(2.*Ay.^2)./D.^3-(2.*Ax.^2)./D.^3).*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(2.*(Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2).*D.^2)+(Ay.*(-(60.*(t(r)-t0).^2)./D.^3+(180.*(t(r)-t0).^3)./D.^4-(120.*(t(r)-t0).^4)./D.^5))./(A_tang.*D.^2)-(2.*Ay.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(A_tang.*D.^3);
        
        Hx(5,2,:) = zeros(size(t));
        Hy(5,2,:) = zeros(size(t));
        Hz(5,2,:) = -(90.*(t(r)-t0).^2)./D.^4+(240.*(t(r)-t0).^3)./D.^5-(150.*(t(r)-t0).^4)./D.^6;
        H(5,2,:) = -(Az.*(-(2.*Az.^2)./D.^3-(2.*Ay.^2)./D.^3-(2.*Ax.^2)./D.^3).*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(2.*(Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2).*D.^2)+(Az.*(-(60.*(t(r)-t0).^2)./D.^3+(180.*(t(r)-t0).^3)./D.^4-(120.*(t(r)-t0).^4)./D.^5))./(A_tang.*D.^2)-(2.*Az.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(A_tang.*D.^3);
        
        
        Hx(1,3,:) = -(60.*(t(r)-t0))./D.^3+(180.*(t(r)-t0).^2)./D.^4-(120.*(t(r)-t0).^3)./D.^5;
        Hy(1,3,:) = zeros(size(t));
        Hz(1,3,:) = zeros(size(t));
        H(1,3,:) = (Ax.*(-(60.*(t(r)-t0))./D.^2+(180.*(t(r)-t0).^2)./D.^3-(120.*(t(r)-t0).^3)./D.^4))./(A_tang.*D.^2);
        
        Hx(2,3,:) = -(90.*(t(r)-t0).^2)./D.^4+(240.*(t(r)-t0).^3)./D.^5-(150.*(t(r)-t0).^4)./D.^6;
        Hy(2,3,:) = zeros(size(t));
        Hz(2,3,:) = zeros(size(t));
        H(2,3,:) = -(Ax.*(-(2.*Az.^2)./D.^3-(2.*Ay.^2)./D.^3-(2.*Ax.^2)./D.^3).*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(2.*(Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2).*D.^2)+(Ax.*(-(60.*(t(r)-t0).^2)./D.^3+(180.*(t(r)-t0).^3)./D.^4-(120.*(t(r)-t0).^4)./D.^5))./(A_tang.*D.^2)-(2.*Ax.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(A_tang.*D.^3);
        
        Hx(3,3,:) = zeros(size(t));
        Hy(3,3,:) = zeros(size(t));
        Hz(2,3,:) = zeros(size(t));
        H(3,3,:) = ((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4)./(A_tang.*D.^2)-(Ax.^2.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./((Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2).*D.^4);
        
        Hx(4,3,:) = zeros(size(t));
        Hy(4,3,:) = zeros(size(t));
        Hz(4,3,:) = zeros(size(t));
        H(4,3,:) = -(Ax.*Ay.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./((Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2).*D.^4);
        
        Hx(5,3,:) = zeros(size(t));
        Hy(5,3,:) = zeros(size(t));
        Hz(5,3,:) = zeros(size(t));
        H(5,3,:) = -(Ax.*Az.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./((Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2).*D.^4);
        
        Hx(1,4,:) = zeros(size(t));
        Hy(1,4,:) = -(60.*(t(r)-t0))./D.^3+(180.*(t(r)-t0).^2)./D.^4-(120.*(t(r)-t0).^3)./D.^5;
        Hz(1,4,:) = zeros(size(t));
        H(1,4,:) = (Ay.*(-(60.*(t(r)-t0))./D.^2+(180.*(t(r)-t0).^2)./D.^3-(120.*(t(r)-t0).^3)./D.^4))./(A_tang.*D.^2);
        
        Hx(2,4,:) = zeros(size(t));
        Hy(2,4,:) = -(90.*(t(r)-t0).^2)./D.^4+(240.*(t(r)-t0).^3)./D.^5-(150.*(t(r)-t0).^4)./D.^6;
        Hz(2,4,:) = zeros(size(t));
        H(2,4,:) = -(Ay.*(-(2.*Az.^2)./D.^3-(2.*Ay.^2)./D.^3-(2.*Ax.^2)./D.^3).*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(2.*(Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2).*D.^2)+(Ay.*(-(60.*(t(r)-t0).^2)./D.^3+(180.*(t(r)-t0).^3)./D.^4-(120.*(t(r)-t0).^4)./D.^5))./(sqrt(Az.^2./D.^2+Ay.^2./D-(Ay.*(-(2.*Az.^2)./D.^3-(2.*Ay.^2)./D.^3-(2.*Ax.^2)./D.^3).*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(2.*(Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2).*D.^2)+(Ay.*(-(60.*(t(r)-t0).^2)./D.^3+(180.*(t(r)-t0).^3)./D.^4-(120.*(t(r)-t0).^4)./D.^5))./(A_tang.*D.^2)-(2.*Ay.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(A_tang.*D.^3).^2+Ax.^2./D.^2).*D.^2)-(2.*Ay.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(A_tang.*D.^3);
        
        Hx(3,4,:) = zeros(size(t));
        Hx(3,4,:) = zeros(size(t));
        Hz(3,4,:) = zeros(size(t));
        H(3,4,:) = -(Ax.*Ay.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./((Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2).*D.^4);
        
        Hx(4,4,:) = zeros(size(t));
        Hy(4,4,:) = zeros(size(t));
        Hz(4,4,:) = zeros(size(t));
        H(4,4,:) = ((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4)./(A_tang.*D.^2)-(Ay.^2.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./((Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2).*D.^4);

        Hx(5,4,:) = zeros(size(t));
        Hy(5,4,:) = zeros(size(t));
        Hz(5,4,:) = zeros(size(t));
        H(5,4,:) = -(Ay.*Az.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./((Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2).*D.^4);
        
        Hx(1,5,:) = zeros(size(t));
        Hy(1,5,:) = zeros(size(t));
        Hz(1,5,:) = -(60.*(t(r)-t0))./D.^3+(180.*(t(r)-t0).^2)./D.^4-(120.*(t(r)-t0).^3)./D.^5;
        H(1,5,:) = (Az.*(-(60.*(t(r)-t0))./D.^2+(180.*(t(r)-t0).^2)./D.^3-(120.*(t(r)-t0).^3)./D.^4))./(A_tang.*D.^2);
        
        Hx(2,5,:) = zeros(size(t));
        Hy(2,5,:) = zeros(size(t));
        Hz(2,5,:) = -(90.*(t(r)-t0).^2)./D.^4+(240.*(t(r)-t0).^3)./D.^5-(150.*(t(r)-t0).^4)./D.^6;
        H(2,5,:) = -(Az.*(-(2.*Az.^2)./D.^3-(2.*Ay.^2)./D.^3-(2.*Ax.^2)./D.^3).*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(2.*(Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2).*D.^2)+(Az.*(-(60.*(t(r)-t0).^2)./D.^3+(180.*(t(r)-t0).^3)./D.^4-(120.*(t(r)-t0).^4)./D.^5))./(sqrt(Az.^2./D.^2+Ay.^2./D-(Az.*(-(2.*Az.^2)./D.^3-(2.*Ay.^2)./D.^3-(2.*Ax.^2)./D.^3).*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(2.*(Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2).*D.^2)+(Az.*(-(60.*(t(r)-t0).^2)./D.^3+(180.*(t(r)-t0).^3)./D.^4-(120.*(t(r)-t0).^4)./D.^5))./(A_tang.*D.^2)-(2.*Az.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(A_tang.*D.^3).^2+Ax.^2./D.^2).*D.^2)-(2.*Az.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./(A_tang.*D.^3);

        Hx(3,5,:) = zeros(size(t));
        Hy(3,5,:) = zeros(size(t));
        Hz(3,5,:) = zeros(size(t));
        H(3,5,:) = -(Ax.*Az.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./((Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2).*D.^4);
        
        Hx(4,5,:) = zeros(size(t));
        Hy(4,5,:) = zeros(size(t));
        Hz(4,5,:) = zeros(size(t));
        H(4,5,:) = -(Ay.*Az.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./((Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2).*D.^4);

        Hx(5,5,:) = zeros(size(t));
        Hy(5,5,:) = zeros(size(t));
        Hz(5,5,:) = zeros(size(t));
        H(5,5,:) = ((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4)./(A_tang.*D.^2)-(Az.^2.*((30.*(t(r)-t0).^2)./D.^2-(60.*(t(r)-t0).^3)./D.^3+(30.*(t(r)-t0).^4)./D.^4))./((Az.^2./D.^2+Ay.^2./D.^2+Ax.^2./D.^2).^(3./2).*D.^4);
    end
end