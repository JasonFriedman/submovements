function [B,J,H] = minimumJerkVelocity1D(t0,D,A,t)
% minimumJerkVelocity1D - evaluate a 1D minimum jerk curve
%
% [B,J,H] = minimumJerkVelocity1D(t0,D,A,t)
%
% see Flash and Hogan (1985) for details on the minimum jerk equation
%
%    t0 = movement start time
%    D  = movement duration
%    Ax = displacement resulting from the movement
%
% The function is evaluated at times t
% Only value of t0 < t < t0 + D will be evaluated (rest will be zero)
%
% The function also optionally returns the first-order and second-order
% partial derivatives, for use with optimization routines
%
% B is the velocity
% J is the gradient (partial derivatives)
% H is the Hessian (second-order partial derivatives)

% Jason Friedman, 2021
% www.curiousjason.com

if numel(t0)>1
    error('t0 must be a number, not a matrix');
end

if numel(D)>1
    error('D must be a number, not a matrix');
end

if numel(A)>1
    error('A must be a number, not a matrix');
end


% nt is normalized time (0 <= nt <= 1)
B = zeros(1,numel(t));
nt = (t-t0)./D;
r = (nt>=0 & nt<=1);
B(r) = A/D * (-60 * nt(r).^3 + 30 * nt(r).^4 + 30 * nt(r).^2);

if nargout > 1
    J = zeros(3,numel(t));
    J(1,r) = A.*(-(60.*(t(r)-t0))./D.^3+(180.*(t(r)-t0).^2)./D.^4-(120.*(t(r)-t0).^3)./D.^5);
    J(2,r) = A.*(-(90.*(t(r)-t0).^2)./D.^4+(240.*(t(r)-t0).^3)./D.^5-(150.*(t(r)-t0).^4)./D.^6);    
    J(3,r) = (30.*(t(r)-t0).^2)./D.^3-(60.*(t(r)-t0).^3)./D.^4+(30.*(t(r)-t0).^4)./D.^5;
        
    if nargout >2
        H(1,1,r) = A.*(60./D.^3-(360.*(t(r)-t0))./D.^4+(360.*(t(r)-t0).^2)./D.^5);      
        H(2,1,r) = A.*((180.*(t(r)-t0))./D.^4-(720.*(t(r)-t0).^2)./D.^5+(600.*(t(r)-t0).^3)./D.^6);
        H(3,1,r) = -(60.*(t(r)-t0))./D.^3+(180.*(t(r)-t0).^2)./D.^4-(120.*(t(r)-t0).^3)./D.^5;
                
        H(1,2,r) = A.*((180.*(t(r)-t0))./D.^4-(720.*(t(r)-t0).^2)./D.^5+(600.*(t(r)-t0).^3)./D.^6);    
        H(2,2,r) = A.*((360.*(t(r)-t0).^2)./D.^5-(1200.*(t(r)-t0).^3)./D.^6+(900.*(t(r)-t0).^4)./D.^7);        
        H(3,2,r) = -(90.*(t(r)-t0).^2)./D.^4+(240.*(t(r)-t0).^3)./D.^5-(150.*(t(r)-t0).^4)./D.^6;
            
        H(1,3,r) = -(60.*(t(r)-t0))./D.^3+(180.*(t(r)-t0).^2)./D.^4-(120.*(t(r)-t0).^3)./D.^5;
        H(2,3,r) = -(90.*(t(r)-t0).^2)./D.^4+(240.*(t(r)-t0).^3)./D.^5-(150.*(t(r)-t0).^4)./D.^6;        
        H(3,3,r) = zeros(size(t(r)));                
    end
end
