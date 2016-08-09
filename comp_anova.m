function [fS,sigS,sig,S] = comp_anova( f, x, K, maxorder, method )
% function to compute anova decomposition (Alg. 1 in [1])
% Author: Steffen Limmer (steffen.limmer@hhi.fraunhofer.de)
% Last update: 28.04.2015
% inputs:
% f: input function, assumed to be integrable using matlab symbolic toolbox
% x: symbolic variables x1,...,xK of function f
% K: number of variables
% maxorder: maximum order of the anova terms, i.e. |S|<=maxorder
% method: method to be speedup symbolic integration
% {'','expand','vpa'}, 'expand' is recommended.
% outputs:
% fS: anova terms, fS(i) depends only on the variables indexed by S(i,:)
% sigS: variances of the anova terms 
% sig: overall variance of the function f
% S: set indexing matrix 

disp('computing anova decomposition for (g \circ f)')

S = fliplr(de2bi(0:2^K-1)); % create binary indexing matrix

idxtmp = find(sum(S,2)<=maxorder); % only compute sets with cardinality <= maxdeg
S = S(idxtmp,:);

for i = 1:size(S,1)
    % symbolic toolbox doesn't support multivariate integration so we use 
    % looping over several single variable integrals
    for j = 1:size(S,2) 
        if j==1
            int_tmp = eval(horzcat(method,'(f)'));
        end       
        if S(i,j) == 0 % index not in u -> integration wrt x(j) must be performed
            int_tmp = int((int_tmp),x(j),0,1,'PrincipalValue',true,'IgnoreSpecialCases',true,'IgnoreAnalyticConstraints',true);
        end
    end
       
    %% identify proper subsets of u
    Stmp = ones(size(S,1),1)*S(i,:) - S;
    idx_psubs = find(sum(Stmp'<0)==0 & sum(Stmp')>=1); % find proper subsets of current index vector U(i,:)
    
    %% compute final anova term
    fS(i) = int_tmp;
    for v=1:numel(idx_psubs) % subtract lower order terms, see example 2.2 (a)
        fS(i) = fS(i) - fS(idx_psubs(v)); %final anova term with index set U(:,i)
    end
    
    %% compute sigma (i.e. power) of anova decomposition terms (this part may be skipped)
    if i == 1      
        %sigu(i) = 0;
    else
        for j = 1:size(S,2)
            % save squared functions for computation of multidimens
            % integral
            if j==1
                int_tmp = eval(horzcat(method,'(fS(i)^2)'));                
            end

            % for variances of decomposition terms (i.e. subset intexed by U(:,i))
            if S(i,j) == 1
                int_tmp = int((int_tmp),x(j),0,1,'PrincipalValue',true,'IgnoreSpecialCases',true,'IgnoreAnalyticConstraints',true);
            else
                int_tmp = int_tmp;
            end     
        end
        
        % variances of complete function (sig) and decomposition terms
        % (sigu)
        sigS(i) = (int_tmp);     
        %sigu(i) = vpa(int_tmp,5);
    end
    
    %% compute sigma of input function
    if i==1               
        for j = 1:size(S,2)
            if j==1
                int_tmp = eval(horzcat(method,'((f-fS(1))^2)'));
            end
            int_tmp = int((int_tmp),x(j),0,1,'PrincipalValue',true,'IgnoreSpecialCases',true,'IgnoreAnalyticConstraints',true);
        end
        sig = (int_tmp);  
    end

disp(horzcat('progress: ',num2str(i/size(S,1))));
end


% resort output according to cardinality of sets
[vals,idx] = sort(sum(S,2));
fS = fS(idx);
sigS = sigS(idx);
S = S(idx,:);

end