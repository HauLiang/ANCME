function x = preliminary_model(A, D, b, K, N)
%
% This code implements the Section IV/A "Preliminary Model" 
%
% The optimization problem may be expressed as
%    minimize   alpha*|| signal - Ax ||_2^2 + || Dx ||_2^2
%
%
% Inputs:
%    A:  dictionary
%    D:  2K block second-order difference matrix
%    signal:  sampled signal
%    K:  number of the modes
%    N:  number of the samples
% Outputs:
%    x:  solution of the above optimization problem
%
% Author: Hao Liang
% Last modified by: 22/12/09
%

% Find a matrix D_orth whose rows are orthogonal to those in D
D_orth = zeros(4*K,2*K*N); ii = 1;
for i = 1:2*K
    D_orth(ii,i*N-1) = 1;
    D_orth(ii+1,i*N) = 1;
    ii = ii+2;
end

% Construct a full rank matrix D_tilde
D_tilde = [D; D_orth];

% Parameter setting
Lambda = A/D_tilde;
Lambda1 = Lambda(:,1:2*K*(N-2));
Lambda2 = Lambda(:,2*K*(N-2)+1:end);
z_matrix = (Lambda2'*Lambda2)\Lambda2';
w_matrix = eye(N,N) - Lambda2*z_matrix;
y = w_matrix*b;
Phi = w_matrix*Lambda1;

% Update w using hierarchical Laplace priors
w = hierarchical_model(Phi,y);

% Update z
z = z_matrix*(b-Lambda1*w);

% Reconstruct x
x = D_tilde\[w;z];

end