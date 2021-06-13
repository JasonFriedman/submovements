function v = lgnbVelocity1D(t0,D,A,sigma,mu,t)
% lgnbVelocity1D - evaluate a 1D lognormal

t1 = t0+D;

% equation from Rohrer & Hogan 2006 Avoiding Spurious ...
v = zeros(1,numel(t));
r = (t>t0 & t<t1);

v(r) = (A * (t1-t0)) ./ (sigma * sqrt(2*pi) * (t(r)-t0) .* (t1 - t(r))) .* ...
    exp( (-1/(2*sigma^2)) * (log((t(r)-t0)./(t1-t(r))) - mu).^2);