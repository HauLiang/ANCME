function w = hierarchical_model(PHI, y)
%
% This code implements the Section IV.B "Hierarchical Model" 
%
% The optimization problem may be expressed as
%    minimize   alpha*|| y - Phi w ||_2^2 + || w ||_2^2
%
% This code is based on the Fast_RVM code available from http://www.miketipping.com/sparsebayes.htm
% [1] Tipping, Michael E, Sparse Bayesian learning and the relevance vector machine, Journal of machine learning research, 2001.
% [2] Tipping, Michael E and Faul, Anita C, Fast marginal likelihood maximisation for sparse Bayesian models, International workshop on artificial intelligence and statistics, 2003
% 
% This code is also based on the paper:
% [3] Babacan S D, Molina R, Katsaggelos A K. Bayesian compressive sensing using Laplace priors. IEEE Transactions on image processing, 2009.
% Please check the accompanying license and the license of [1] and [2] before using. 
%
% Inputs:
%    Phi:  transformed dictionary
%    y:  transformed sampled signal
% Outputs:
%    x:  solution of the above optimization problem
%
% Author: Hao Liang
% Last modified by: 22/12/09
%

% Parameter setting
eta = 1e-8;  % threshold for stopping
inv_beta =  std(y)^2/1e2;  % initial noise variance
maxIter = 1000;
deleted = [];
 
% Find initial alpha
[~,M] = size(PHI);
PHIy = PHI'*y;
PHI2 = sum(PHI.^2)';
ratio = (PHIy.^2)./PHI2;
[maxr,index] = max(ratio);
inv_alpha = PHI2(index)/(maxr-inv_beta);

% Compute initial mu, Sig, S, Q
phi = PHI(:,index);
Hessian = inv_alpha + phi'*phi/inv_beta;
Sig = 1/Hessian;
mu = Sig*PHIy(index)/inv_beta;
left = PHI'*phi/inv_beta;
S = PHI2/inv_beta-Sig*left.^2;
Q = PHIy/inv_beta-Sig*PHIy(index)/inv_beta*left;

for count = 1:maxIter
    
    % Calculate si and qi
    s = S; q = Q;
    s(index) = inv_alpha.*S(index)./(inv_alpha-S(index));
    q(index) = inv_alpha.*Q(index)./(inv_alpha-S(index));
    
    % Eq.(43)
    lambda = (2*length(index) - 2) / (sum(1./inv_alpha));

    A = lambda + s - q.^2;
    B = 2*lambda.*s + s.^2;
    C = lambda.*s.^2;
        
    theta = q.^2-s;
    discriminant = B.^2 - 4.*A.*C;
    nextAlphas = (-B - sqrt(discriminant) ) ./ (2*A);

    % Choose the next alpha that maximizes marginal likelihood
    ml = -inf*ones(1,M);
    ig0 = find(theta>lambda);
    
    % Indices for reestimation
    [ire,~,which] = intersect(ig0,index);
    if ~isempty(ire)
        Alpha = nextAlphas(ire);
        delta = (inv_alpha(which)-Alpha)./(Alpha.*inv_alpha(which));
        ml(ire) = q(ire).^2./ (Alpha + s(ire)) + log(Alpha ./ (Alpha + s(ire))) - lambda./ Alpha ...
            -q(ire).^2./ (inv_alpha(which) + s(ire)) - log(inv_alpha(which) ./ (inv_alpha(which) + s(ire))) + lambda./ inv_alpha(which);
    end
    
    % Indices for adding
    iad = setdiff(ig0,ire);
    if ~isempty(iad)
        Alpha = nextAlphas(iad);
        ml(iad) = log(Alpha ./ (Alpha + s(iad)) )+ q(iad).^2 ./ (Alpha + s(iad)) - lambda./Alpha;
        which = intersect(deleted,iad);
        ml(which) = -inf;   
    end
    is0 = setdiff((1:M),ig0);
    
    % Indices for deleting
    [ide,~,which] = intersect(is0,index);
    if ~isempty(ide)
         if length(index) == 1
             ml(ide) = -inf;
         else
             ml(ide) = -q(ide).^2 ./ (inv_alpha(which) + s(ide)) - log( inv_alpha(which) ./(inv_alpha(which) + s(ide))) + lambda./inv_alpha(which);
         end
    end

     % check convergence
     [ML(count),idx] = max(ml);
    if count > 2 && abs(ML(count)-ML(count-1)) < abs(ML(count)-ML(1))*eta
        break;
    end

    % Update alphas
    which = find(index==idx);
    if theta(idx) > lambda
        if ~isempty(which)  % re-estimate a basis
          
            Alpha = nextAlphas(idx);
            Sigii = Sig(which,which); mui = mu(which); Sigi = Sig(:,which);
            delta = Alpha-inv_alpha(which);
            ki = delta/(1+Sigii*delta);
            mu = mu-ki*mui*Sigi;
            Sig = Sig-ki*(Sigi*Sigi');
            comm = PHI'*(phi*Sigi)/inv_beta;
            S = S + ki*comm.^2;
            Q = Q + ki*mui*comm;
            inv_alpha(which) = Alpha;
            
        else     % add a basis
            
            Alpha = nextAlphas(idx);
            phii = PHI(:,idx); Sigii = 1/(Alpha+S(idx)); mui = Sigii*Q(idx);
            comm1 = Sig*(phi'*phii)/inv_beta;
            ei = phii-phi*comm1;
            off = -Sigii*comm1;
            Sig = [Sig+Sigii*(comm1*comm1'), off; off', Sigii];
            mu = [mu-mui*comm1; mui];
            comm2 = PHI'*ei/inv_beta;
            S = S - Sigii*comm2.^2;
            Q = Q - mui*comm2;
            
            index = [index;idx];
            inv_alpha = [inv_alpha;Alpha];
            phi = [phi,phii];
            
        end
    else
        if ~isempty(which) && length(index) > 1  % delete a basis
            
            deleted = [deleted idx];
            Sigii = Sig(which,which); mui = mu(which); Sigi = Sig(:,which);
            Sig = Sig-Sigi*Sigi'/Sigii; Sig(:,which) = []; Sig(which,:) = [];
            mu  = mu-mui/Sigii*Sigi; mu(which) = [];
            comm = PHI'*(phi*Sigi)/inv_beta;
            S = S + comm.^2/Sigii;
            Q = Q + mui/Sigii*comm;
            
            index(which) = [];
            inv_alpha(which) = [];
            phi(:,which) = [];
        elseif ~isempty(which) && length(index) == 1
            % Something is wrong, trying to delete the only coefficient that has been added.
            break;
        end
            
    end    
end

% Estimation of w
weights	= mu; used = index;
w = zeros(M,1);  w(used) = weights;

end