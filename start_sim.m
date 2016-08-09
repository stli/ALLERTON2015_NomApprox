% Script to generate figures 2(a)-(f)\(b) in [1]
% Author: Steffen Limmer (steffen.limmer@hhi.fraunhofer.de)
% Last update: 17.09.2015
% References: 
% [1] @article{Li_2015,
%	author = "Limmer, S. and Mohammadi, J. and Stanczak, S.",
%	title = "A Simple Algorithm for Nomographic Approximation",
%	year = "2015"}

%% setup simulation parameters
K = 2;  % number of variables

%% define function and symbolic variables
% IMPORTANT2: SYMBOLIC TOOLBOX IS REQUIRED
x = sym('x', [K,1]);    % define symbolic variables
x = sym(x, 'positive'); % assume x is positive and domain is [0,1]^d
% paper example
f = 1/9*(x(1)+x(2)+x(1)*x(2))^2; % f:[0,1]^K \to [0,1]

D = 20; % degree of polynomial psi^{-1}

%% define polynomial p:=psi^{-1}
t = sym('t');           % variable for polynomial
g = 0;                  % initialize p
z = sym('z', [D 1]);    % initialize a 
for d = 1:D
    g = g + z(d)*t^d;   % polynomial of degree D
end

%% compute matrices A,B,M
tic
[A,B,M] = gen_matrices((f),x,K,D,'expand');
toc

%% compute optimal coefficients of p using SDR
P = inv((M));
delta = 1e-3;
[zopt, sdpval, rlq] = opt_sdr(vpa(A),vpa(B),vpa(M),D,delta,vpa(P));

%% compute anova decomposition of phi := p \circ f 
popt = subs(g,z,double(zopt)); 
maxdeg = 1; % maximum degree for the anova terms
[phiS,sigS,sig,S] = comp_anova( subs(popt,t,f),x,K,maxdeg,'expand' );

%% compute numerical inverse psi := p^{-1}(\xi + phi0) with dom(psi)=[0,1] 
psi = comp_numinverse(popt,double(phiS(1))); % \psi:\Omega' \to [0,1]

%% nomographic approximation is obtained by 
% fhat = psi( sum( phiS(find(sum(S,2)<=1)) );

%% plot results 2(a)-(f)\(b)
% plot with larger domain to illustrate inner approximation effects
% if K=2 plot resulting functions and errors
psi_minin = double( subs(popt,0)-phiS(1) )-0.02;
psi_maxin = double( subs(popt,1)-phiS(1) )+0.02;

plot_results(f, x, psi, [0,phiS(2:end)], S, maxdeg, K, psi_minin, psi_maxin);