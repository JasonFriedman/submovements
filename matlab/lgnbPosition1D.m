function x = lgnbPosition1D(t0,D,A,sigma,mu,x0,t)
% lgnbPosition1D - evaluate a 1D lognormal position

% find by integration
v = lgnbVelocity1D(t0,D,A,sigma,mu,t);
x = cumtrapz(t,v) + x0;