function [ ] = plot_results( f, x, psi, phiS, S, maxdeg, K, ymin, ymax )
% function to plot fiugres 2(a)-(f)\(b) in [1]
% Author: Steffen Limmer (steffen.limmer@hhi.fraunhofer.de)
% Last update: 17.09.2015
% References: 
% [1] @article{Li_2015,
%	author = "Limmer, S. and Mohammadi, J. and Stanczak, S.",
%	title = "A Simple Algorithm for Nomographic Approximation",
%	year = "2015"}

f_matfun = matlabFunction(f,'vars',{x}); % speedup using matlabfunction

%% define functions phi_emptyset, phi_1, phi_2, only applies to plotting of 2D functions
if K==2
    phi1_matfun = matlabFunction( vpa( phiS( find(ismember(S,[1,0],'rows')) ) ),'vars',{x});
    phi2_matfun = matlabFunction( vpa( phiS( find(ismember(S,[0,1],'rows')) ) ),'vars',{x});
end

phihat_sym = vpa( sum( phiS(find(sum(S,2)<=maxdeg)) ) );
phihat_matfun = matlabFunction(phihat_sym,'vars',{x});

%% plot psi
y = linspace(ymin,ymax,1e3);
for l=1:numel(y)
    yi(l) = psi(y(l));
end
figure
plot(y,yi)
xlabel('y')
ylabel('psi(y)')
title('numerical inverse psi')


%% compute approximation error
if K~=2
    disp('note: dimensions do not match. skip plotting of results...')
else
    [xn1,xn2] = meshgrid(0:0.02:1);
    ngrid = 0:0.02:1;
    for i = 1:size(xn1,1)
        for j = 1:size(xn1,1)
            xn = [xn1(i,j);xn2(i,j)];
            fref(i,j) = f_matfun(xn);

            %only applies to plotting of 2D functions
            phi1(i,j) = phi1_matfun(xn);
            phi2(i,j) = phi2_matfun(xn);


            % nomographic approximation
            phihat(i,j) = phihat_matfun(xn);
            fhat(i,j) = double( psi( phihat(i,j) ) );
            e(i,j) = (fref(i,j) - fhat(i,j));
        end
    end

    %% plot reference function
    figure
    mesh(ngrid,ngrid,fref)
    xlabel('x1')
    ylabel('x2')
    zlabel('zlabel')
    title('f evaluated for [0,1]^2')

    %% plot approximation error
    figure
    mesh(ngrid,ngrid,fhat)
    xlabel('x1')
    ylabel('x2')
    zlabel('zlabel')
    title('nomographic approximation')
    
    %% plot approximation error
    figure
    mesh(ngrid,ngrid,e)
    xlabel('x1')
    ylabel('x2')
    zlabel('zlabel')
    title('nomographic approximation error')

    %% plot \phi_1
    figure
    mesh(ngrid,ngrid,phi1)
    xlabel('x1')
    ylabel('x2')
    zlabel('zlabel')
    title('\phi_1 evaluated for [0,1]^2')

    %% plot \phi_2
    figure
    mesh(ngrid,ngrid,phi2)
    xlabel('x1')
    ylabel('x2')
    zlabel('zlabel')
    title('\phi_2 evaluated for [0,1]^2')
    end

end