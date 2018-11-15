function [Q, R] = mgs(A)
% MGS Direct implementation of Algorithm 8.1
% Input: an mxn matrix A
% Outputs: matrix Q with orthonormal columns, upper triangular matrix R

    % Extract size of a and create output matrices
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);
    
    % Assignment 
    V = A; % create the v vectors
    
    % Main loop
    for i=1:n
        R(i,i) = norm(V(:,i)); % calculate r_(ii)
        Q(:, i) = V(:,i)/R(i,i); % create q_i
        for j=i+1:n
           R(i,j) = Q(:,i)'*V(:,j); % projection step
           V(:,j) = V(:,j) - R(i,j)*Q(:,i); % updating
        end
    end

end