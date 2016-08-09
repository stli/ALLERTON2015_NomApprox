function [ A,B,M ] = gen_matrices( f,x,K,D,method )
% function to generate matrices A,B,M in [1]
% Author: Steffen Limmer (steffen.limmer@hhi.fraunhofer.de)
% Last update: 17.09.2015
% inputs:
% f : input function
% x : symbolic variables x1,...,xK of function f
% K: number of variables
% D: degree of p -> A,B,M \in R^{D \times D}
% method: method to be speedup symbolic integration
% {'','expand','vpa'}, 'expand' is recommended.
% outputs:
% A: a'*A*a = sum(sigma_k(p(a) \circ f)), obtained by eq (13-14)
% B: a'*B*a = sigma(p(a) \circ f), obtained by eq (12)
% M: change of basis matrix, obtained by eq (11)

disp('computing matrices A B M...')
%% compute matrix M as described in Lemma 2 and (11) in [1]
for k = 0:D-1
    for r = 0:k
        Mt(k+1,r+1) = sym( nchoosek(k,r)/nchoosek(D-1,r) );
    end
end
M = Mt * diag([1:D]);

%% compute matrix B as described in (12)
b = sym(0); % enforce entries of b to be symbolic
for d = 1:2*D
    int_tmp = eval(horzcat(method,'(f^d)'));
    for k=1:K
        int_tmp = ( int(int_tmp,x(k),0,1,'IgnoreSpecialCases',true,'IgnoreAnalyticConstraints',true) );
    end
    b(d) = int_tmp;
    disp(horzcat('progress: ',num2str(d/(2*D))));
end
B1 = hankel(b(2:D+1),b(D+1:2*D));
B2 = b(1:D)'*b(1:D);
B = B1 - B2;

%% compute matrix A as described in (13),(14)
A1k = sym(zeros(D,D,K)); % enforce entries of Ak to be symmetric
A = sym(zeros(D,D)); % enforce entries of A to be symmetric
for k = 1:K
    for i = 1:D
        for j = 1:i % matrices are symmetric so we only compute lower triangular part
            int_tmp1 = (f^i);
            int_tmp2 = (f^j);
            for kk = setdiff(1:K,k) % compute integrals w.r.t. variables {1,...K}\k
                int_tmp1 = ( int(int_tmp1,x(kk),0,1,'IgnoreSpecialCases',true,'IgnoreAnalyticConstraints',true) );
                int_tmp2 = ( int(int_tmp2,x(kk),0,1,'IgnoreSpecialCases',true,'IgnoreAnalyticConstraints',true) );
            end
            A1k(i,j,k) = int( eval(horzcat(method,'(int_tmp1*int_tmp2)')), x(k),0,1,'PrincipalValue',true,'IgnoreSpecialCases',true,'IgnoreAnalyticConstraints',true);
        end
    end
    A1k(:,:,k) = A1k(:,:,k) + A1k(:,:,k)' - diag(diag(A1k(:,:,k))); % compute A^{(1)}(k)
    disp(horzcat('progress: ',num2str(k/(K))));  
    A = A + A1k(:,:,k) - B2; % compute final matrix A = sum( Ak - B2 )
end

end

