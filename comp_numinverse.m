function [ psi ] = comp_numinverse( popt, phi0 )
% function to compute numerical inverse
% Author: Steffen Limmer (steffen.limmer@hhi.fraunhofer.de)
% Last update: 17.09.2015
% inputs: 
% popt: optimal polynomial (:= psi^{-1}):= polynomial in matlab symbolic form
% output: 
% psi: inverse of p, dom(psi) = range(p); range(psi) = [0,1]

xvals = linspace(0,1,1e4); % range(f) = [0,1];
pmatfun = matlabFunction(vpa(popt)); % generate matlab function for faster function calls
pvals = pmatfun( xvals );

%% projection of y onto [min(range(p)),max(range(p))], values outside this
% interval are caused by truncation of anova decomposition
% matlab function for faster evaluation
proj = @(y,a,b) min( max(y,a), b);

%% find minimum between x and xvec, output is index in xvec with minimum distance
findmin = @(x,xvec) find((abs(x-xvec) - min(abs(x-xvec)) == 0 ),1); 

%% function to compute psi:=p^{-1} using functions defined above
psi = @(y) xvals( findmin(proj(y+phi0,pvals(1),pvals(end)),pvals) );

end
