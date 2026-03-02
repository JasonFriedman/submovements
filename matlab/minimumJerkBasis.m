function [nt,r,before,after,positionBasis,velocityBasis] = minimumJerkBasis(t0,D,t)
% MINIMUMJERKBASIS - shared normalized-time basis for minimum jerk functions

nt = (t-t0)./D;
r = (nt>=0 & nt<=1);
before = nt<0;
after = nt>1;
positionBasis = (-15 * nt.^4 + 6 * nt.^5 + 10 * nt.^3);
velocityBasis = (-60 * nt.^3 + 30 * nt.^4 + 30 * nt.^2);
