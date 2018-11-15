function [W, R] = house(A)
% HOUSE Compute an implicit representation of a
% full QR factorization for an mxn matrix A using 
% Householder reflections
% Inputs: mxn matrix A, m>=n
% Ouputs: lower triangular matrix W housing vectors that define
%         the reflections; upper triangular matrix R

    % Extract size of A and create output matrices
    [m, n] = size(A);
    W = zeros(m, n);
    R = A; % initialize R to A
    
    % Main loop
    for k=1:n
        x= R(k:m,k);
        signX1 = sign(x(1)); % define the sign
        % Use the normalization from the lecture for sign of a zero first
        % element
        if signX1==0
            signX1= 1;
        end
        I = eye(m-k+1); % create identity to obtain e1
        W(k:m,k) = signX1*norm(x)*I(:,1) + x; % create v_k
        W(:,k) = W(:,k)/norm(W(:,k)); % normalize v_k
        R(k:m,k:n) = R(k:m,k:n)- 2*W(k:m,k)*( W(k:m,k)'*R(k:m,k:n)); % update R
    end
end