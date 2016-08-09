function [ zopt, sdpval, rlq ] = opt_sdr( A,B,M,D,delta,P )
% function to compute zopt using semidefinite relaxation of eq (16) in [1]
% Author: Steffen Limmer (steffen.limmer@hhi.fraunhofer.de)
% Last update: 28.04.2015
% zopt is solution to min(z'*A*z)/(z'*B*z) s.t. M*z>=zeros(D,1)
% z~=zeros(D,1)
% Author: Steffen Limmer (steffen.limmer@hhi.fraunhofer.de)
% Last update: 28.04.2015
% inputs:
% A \in R^{D \times D}, matrix in the numerator, obtained by eq (14-15)
% B \in R^{D \times D}, matrix in the denominator, obtained by eq (13)
% M \in R^{D \times D}, change of basis matrix, obtained by eq (11-12)
% D: size of zopt
% delta: bound for trace(B*Z)
% P: preconditioning matrix
% outputs:
% zopt: solution to optimization problem obtained by semidefinite
% relaxation

% check condition number of involved matrices
disp('solving semidefinite optimization problem...');

%% solve semidefinite relaxation of (18)
% requires cvx (http://cvxr.com/), testet with ver 2.1 (Build 1085)
if nargin < 6
    if rcond(A) <= 1e-10 || rcond(B) <= 1e-10
        disp('matrices are ill-conditioned, using matrix preconditioning');
        P = inv(double(M) ); %preconditioning
    else
        P = eye(D);
    end
elseif nargin < 5
    delta = 1e-2;
end

cvx_begin sdp quiet
    cvx_solver sedumi
    cvx_precision best

    variable Z(D,D) symmetric; 
    
    maximize ( trace( double( P'*A*P )*Z ) )

    subject to
    trace( double( P'*B*P )*Z ) == delta
    Z == semidefinite(D) 
    for i=1:D
        for j=1:i
            trace( Z*double( P'*M(j,:)'*M(i,:)*P ) ) >=0
        end
    end 
cvx_end


% try to recover zopt by rank 1 approximation
sdpval = cvx_optval/delta;  % objective value of the sdp
[Q S] = eigs(double(Z));            % obtain eigenvalue decomposition for rank 1 approximation
zsdp = P*Q(:,1);

if abs(S(2,2)/S(1,1))<=1e-9 && all( M*zsdp >= 0 ) % [sign(M*zsdp)]_1=sign(zsdp(1))
    disp('rank and cone constraint fulfilled, optimal solution found');
    zopt = zsdp;
elseif abs(S(2,2)/S(1,1))<=1e-9 && all( M*zsdp <= 0 ) % [sign(M*zsdp)]_1=sign(zsdp(1))
    disp('rank and cone constraint fulfilled, optimal solution found');
    zopt = -zsdp;
else
    disp('warning: cone or rank constraint not fulfilled. projecting onto constraints');
    if all( M*zsdp >= 0 ) 
        zopt = zsdp; % pick positive solution
    elseif all( M*zsdp <= 0)
        zopt = -zsdp; % pick negative solution
    else % project onto cone cone K := {z | M*z >= 0}
        rlqf = @(z) double(z'*A*z)/double(z'*B*z);
        zopt_tmp(:,1) = double(inv(M)*pos( double(M*zsdp) ) );
        zopt_tmp(:,2) = double(inv(M)*pos( double(M*(-zsdp)) ) );
        [~,idx] = max([rlqf(zopt_tmp(:,1)),rlqf(zopt_tmp(:,2))]);
        zopt = zopt_tmp(:,idx );
    end
end

if any( ~isfinite( zopt ) )==1 %%check if zopt is valid
    warning('error: zopt is not finite');
end

% compute objective value for rank 1 approx
% needs inversion of diagonal scaling matrix Dd
rlq = double(zopt'*A*zopt)/double(zopt'*B*zopt);

% check tightness if sdp is tight
disp(horzcat('semidefinite program objective value: ',num2str(sdpval)));
disp(horzcat('rayleigh quotient objective value: ',num2str(rlq)));

if rlq < sdpval-1e-2
    disp('warning: semidefinite relaxation is not tight');
end

end