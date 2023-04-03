function [Bx,By,B,Jx,Jy,J,Hx,Hy,H] = minimumJerkVelocity2D(t0,D,Ax,Ay,t)
% minimumJerkVelocity2D - evaluate a minimum jerk velocity curve with seperate displacement for x / y
%
% see Flash and Hogan (1985) for details on the minimum jerk equation
%
%    t0 = movement start time
%    D  = movement duration
%    Ax = displacement resulting from the movement (x)
%    Ay = displacement resulting from the movement (y)
%
% The function is evaluated at times t
%
% The function also optionally returns the first-order and second-order
% partial derivatives, for use with optimization routines
%
% Bx, By and B are the x velocity, y velocity and tangential velocities
% Jx, Jy and J are the gradients (partial derivatives) of the same quantities
% Hx, Hy and H are the Hessian (second-order partial derivatives)
%
% [Bx,By,B,Jx,Jy,J,Hx,Hy,H] = minimumJerkVelocity2D(t0,D,Ax,Ay,t)

% Jason Friedman, 2023
% www.curiousjason.com

if numel(t0)>1
    error('t0 must be a number, not a matrix');
end

if numel(D)>1
    error('D must be a number, not a matrix');
end

if numel(Ax)>1
    error('Ax must be a number, not a matrix');
end

if numel(Ay)>1
    error('Ay must be a number, not a matrix');
end


% nt is normalized time (0 <= nt <= 1)
Bx = zeros(1,numel(t));
By = zeros(1,numel(t));
B = zeros(1,numel(t));
nt = (t-t0)./D;
r = (nt>=0 & nt<=1);

Bx(r) = Ax/D * (-60 * nt(r).^3 + 30 * nt(r).^4 + 30 * nt(r).^2);
By(r) = Ay/D * (-60 * nt(r).^3 + 30 * nt(r).^4 + 30 * nt(r).^2);

A_tang = sqrt((Ax/D).^2 + (Ay/D).^2);
B(r) = A_tang * (-60 * nt(r).^3 + 30 * nt(r).^4 + 30 * nt(r).^2);

if nargout > 3
    Jx = zeros(4,numel(t));
    Jy = zeros(4,numel(t));
    J = zeros(4,numel(t));

    % These are calculated by derivegradianthessian.m
    Jx(1,r) = -Ax.*(1.0./D.^3.*(t(r).*2.0-t0.*2.0).*3.0e+1-1.0./D.^4.*(t(r)-t0).^2.*1.8e+2+1.0./D.^5.*(t(r)-t0).^3.*1.2e+2);
    Jx(2,r) = -Ax.*(1.0./D.^4.*(t(r)-t0).^2.*9.0e+1-1.0./D.^5.*(t(r)-t0).^3.*2.4e+2+1.0./D.^6.*(t(r)-t0).^4.*1.5e+2);
    Jx(3,r) = 1.0./D.^3.*(t(r)-t0).^2.*3.0e+1-1.0./D.^4.*(t(r)-t0).^3.*6.0e+1+1.0./D.^5.*(t(r)-t0).^4.*3.0e+1;

    Jy(1,r) = -Ay.*(1.0./D.^3.*(t(r).*2.0-t0.*2.0).*3.0e+1-1.0./D.^4.*(t(r)-t0).^2.*1.8e+2+1.0./D.^5.*(t(r)-t0).^3.*1.2e+2);
    Jy(2,r) = -Ay.*(1.0./D.^4.*(t(r)-t0).^2.*9.0e+1-1.0./D.^5.*(t(r)-t0).^3.*2.4e+2+1.0./D.^6.*(t(r)-t0).^4.*1.5e+2);

    Jy(4,r) = 1.0./D.^3.*(t(r)-t0).^2.*3.0e+1-1.0./D.^4.*(t(r)-t0).^3.*6.0e+1+1.0./D.^5.*(t(r)-t0).^4.*3.0e+1;
    
    J(1,r) = (1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^3.*(D-t(r)+t0).^4.*4.0-1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^3.*4.0).*1.0./sqrt(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).*-1.5e+1;
    J(2,r) = (1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^3.*4.0-1.0./D.^11.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4.*1.0e+1).*1.0./sqrt(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).*1.5e+1;
    J(3,r) = Ax.*1.0./D.^10.*(t(r)-t0).^4.*(D-t(r)+t0).^4.*1.0./sqrt(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).*3.0e+1;
    J(4,r) = Ay.*1.0./D.^10.*(t(r)-t0).^4.*(D-t(r)+t0).^4.*1.0./sqrt(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).*3.0e+1;

    if nargout >6
        Hx = zeros(4,4,numel(t));
        Hy = zeros(4,4,numel(t));
        H = zeros(4,4,numel(t));
        
        Hx(1,1,r) = Ax.*(1.0./D.^4.*(t(r).*2.0-t0.*2.0).*-1.8e+2+1.0./D.^5.*(t(r)-t0).^2.*3.6e+2+1.0./D.^3.*6.0e+1);
        Hx(1,2,r) = Ax.*(1.0./D.^4.*(t(r).*2.0-t0.*2.0).*9.0e+1-1.0./D.^5.*(t(r)-t0).^2.*7.2e+2+1.0./D.^6.*(t(r)-t0).^3.*6.0e+2);
        Hx(1,3,r) = 1.0./D.^3.*(t(r).*2.0-t0.*2.0).*-3.0e+1+1.0./D.^4.*(t(r)-t0).^2.*1.8e+2-1.0./D.^5.*(t(r)-t0).^3.*1.2e+2;

        Hx(2,1,r) = Ax.*(1.0./D.^4.*(t(r).*2.0-t0.*2.0).*9.0e+1-1.0./D.^5.*(t(r)-t0).^2.*7.2e+2+1.0./D.^6.*(t(r)-t0).^3.*6.0e+2);
        Hx(2,2,r) = Ax.*(1.0./D.^5.*(t(r)-t0).^2.*3.6e+2-1.0./D.^6.*(t(r)-t0).^3.*1.2e+3+1.0./D.^7.*(t(r)-t0).^4.*9.0e+2);
        Hx(2,3,r) = 1.0./D.^4.*(t(r)-t0).^2.*-9.0e+1+1.0./D.^5.*(t(r)-t0).^3.*2.4e+2-1.0./D.^6.*(t(r)-t0).^4.*1.5e+2;

        Hx(3,1,r) = 1.0./D.^3.*(t(r).*2.0-t0.*2.0).*-3.0e+1+1.0./D.^4.*(t(r)-t0).^2.*1.8e+2-1.0./D.^5.*(t(r)-t0).^3.*1.2e+2;
        Hx(3,2,r) = 1.0./D.^4.*(t(r)-t0).^2.*-9.0e+1+1.0./D.^5.*(t(r)-t0).^3.*2.4e+2-1.0./D.^6.*(t(r)-t0).^4.*1.5e+2;

        Hy(1,1,r) = Ay.*(1.0./D.^4.*(t(r).*2.0-t0.*2.0).*-1.8e+2+1.0./D.^5.*(t(r)-t0).^2.*3.6e+2+1.0./D.^3.*6.0e+1);
        Hy(1,2,r) = Ay.*(1.0./D.^4.*(t(r).*2.0-t0.*2.0).*9.0e+1-1.0./D.^5.*(t(r)-t0).^2.*7.2e+2+1.0./D.^6.*(t(r)-t0).^3.*6.0e+2);

        Hy(1,4,r) = 1.0./D.^3.*(t(r).*2.0-t0.*2.0).*-3.0e+1+1.0./D.^4.*(t(r)-t0).^2.*1.8e+2-1.0./D.^5.*(t(r)-t0).^3.*1.2e+2;
        Hy(2,1,r) = Ay.*(1.0./D.^4.*(t(r).*2.0-t0.*2.0).*9.0e+1-1.0./D.^5.*(t(r)-t0).^2.*7.2e+2+1.0./D.^6.*(t(r)-t0).^3.*6.0e+2);
        Hy(2,2,r) = Ay.*(1.0./D.^5.*(t(r)-t0).^2.*3.6e+2-1.0./D.^6.*(t(r)-t0).^3.*1.2e+3+1.0./D.^7.*(t(r)-t0).^4.*9.0e+2);

        Hy(2,4,r) = 1.0./D.^4.*(t(r)-t0).^2.*-9.0e+1+1.0./D.^5.*(t(r)-t0).^3.*2.4e+2-1.0./D.^6.*(t(r)-t0).^4.*1.5e+2;

        Hy(4,1,r) = 1.0./D.^3.*(t(r).*2.0-t0.*2.0).*-3.0e+1+1.0./D.^4.*(t(r)-t0).^2.*1.8e+2-1.0./D.^5.*(t(r)-t0).^3.*1.2e+2;
        Hy(4,2,r) = 1.0./D.^4.*(t(r)-t0).^2.*-9.0e+1+1.0./D.^5.*(t(r)-t0).^3.*2.4e+2-1.0./D.^6.*(t(r)-t0).^4.*1.5e+2;

        H(1,1,r) = 1.0./D.^20.*(Ax.^2+Ay.^2).^2.*(t(r)-t0).^6.*(D-t(r)+t0).^6.*1.0./(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).^(3.0./2.0).*(D.*t(r).*-6.0+D.*t0.*6.0-t(r).*t0.*1.2e+1+D.^2+t(r).^2.*6.0+t0.^2.*6.0).*6.0e+1;
        H(1,2,r) = 1.0./D.^21.*(Ax.^2+Ay.^2).^2.*(t(r)-t0).^7.*(D-t(r)+t0).^6.*1.0./(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).^(3.0./2.0).*(D.*t(r).*-1.2e+1+D.*t0.*1.2e+1-t(r).*t0.*2.0e+1+D.^2.*3.0+t(r).^2.*1.0e+1+t0.^2.*1.0e+1).*6.0e+1;
        H(1,3,r) = Ax.*1.0./D.^20.*(Ax.^2+Ay.^2).*(t(r)-t0).^7.*(D-t(r).*2.0+t0.*2.0).*(D-t(r)+t0).^7.*1.0./(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).^(3.0./2.0).*-6.0e+1;
        H(1,4,r) = Ay.*1.0./D.^20.*(Ax.^2+Ay.^2).*(t(r)-t0).^7.*(D-t(r).*2.0+t0.*2.0).*(D-t(r)+t0).^7.*1.0./(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).^(3.0./2.0).*-6.0e+1;
        H(2,1,r) = 1.0./D.^21.*(Ax.^2+Ay.^2).^2.*(t(r)-t0).^7.*(D-t(r)+t0).^6.*1.0./(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).^(3.0./2.0).*(D.*t(r).*-1.2e+1+D.*t0.*1.2e+1-t(r).*t0.*2.0e+1+D.^2.*3.0+t(r).^2.*1.0e+1+t0.^2.*1.0e+1).*6.0e+1;
        H(2,2,r) = 1.0./D.^22.*(Ax.^2+Ay.^2).^2.*(t(r)-t0).^8.*(D-t(r)+t0).^6.*1.0./(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).^(3.0./2.0).*(D.*t(r).*-2.0e+1+D.*t0.*2.0e+1-t(r).*t0.*3.0e+1+D.^2.*6.0+t(r).^2.*1.5e+1+t0.^2.*1.5e+1).*6.0e+1;
        H(2,3,r) = Ax.*1.0./D.^21.*(Ax.^2+Ay.^2).*(t(r)-t0).^8.*(D.*3.0-t(r).*5.0+t0.*5.0).*(D-t(r)+t0).^7.*1.0./(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).^(3.0./2.0).*-3.0e+1;
        H(2,4,r) = Ay.*1.0./D.^21.*(Ax.^2+Ay.^2).*(t(r)-t0).^8.*(D.*3.0-t(r).*5.0+t0.*5.0).*(D-t(r)+t0).^7.*1.0./(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).^(3.0./2.0).*-3.0e+1;
        H(3,1,r) = Ax.*1.0./D.^20.*(Ax.^2+Ay.^2).*(t(r)-t0).^7.*(D-t(r).*2.0+t0.*2.0).*(D-t(r)+t0).^7.*1.0./(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).^(3.0./2.0).*-6.0e+1;
        H(3,2,r) = Ax.*1.0./D.^21.*(Ax.^2+Ay.^2).*(t(r)-t0).^8.*(D.*3.0-t(r).*5.0+t0.*5.0).*(D-t(r)+t0).^7.*1.0./(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).^(3.0./2.0).*-3.0e+1;
        H(3,3,r) = Ay.^2.*1.0./D.^20.*(t(r)-t0).^8.*(D-t(r)+t0).^8.*1.0./(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).^(3.0./2.0).*3.0e+1;
        H(3,4,r) = Ax.*Ay.*1.0./D.^20.*(t(r)-t0).^8.*(D-t(r)+t0).^8.*1.0./(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).^(3.0./2.0).*-3.0e+1;
        H(4,1,r) = Ay.*1.0./D.^20.*(Ax.^2+Ay.^2).*(t(r)-t0).^7.*(D-t(r).*2.0+t0.*2.0).*(D-t(r)+t0).^7.*1.0./(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).^(3.0./2.0).*-6.0e+1;
        H(4,2,r) = Ay.*1.0./D.^21.*(Ax.^2+Ay.^2).*(t(r)-t0).^8.*(D.*3.0-t(r).*5.0+t0.*5.0).*(D-t(r)+t0).^7.*1.0./(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).^(3.0./2.0).*-3.0e+1;
        H(4,3,r) = Ax.*Ay.*1.0./D.^20.*(t(r)-t0).^8.*(D-t(r)+t0).^8.*1.0./(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).^(3.0./2.0).*-3.0e+1;
        H(4,4,r) = Ax.^2.*1.0./D.^20.*(t(r)-t0).^8.*(D-t(r)+t0).^8.*1.0./(1.0./D.^10.*(Ax.^2+Ay.^2).*(t(r)-t0).^4.*(D-t(r)+t0).^4).^(3.0./2.0).*3.0e+1;
    end
end
