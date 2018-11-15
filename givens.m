function R = givens(A)
% givens Computes the matrix R of the QR 
% decomposition of a matrix A using Givens rotations
% Inputs: mxn matrix A, m>=n
% Ouputs: upper triangular matrix R  
    
    % Extract size of A and create output matrices
    [m, n] = size(A);
    R = A; % initialize R to A
    
    % Main loop
    for k=1:n % go over columns
        for i=m-1:-1:k
           if norm(R(i:i+1,k))~=0  % check so that we don't divide by zero
                r = norm(R(i:i+1,k));
                s = R(i+1,k)/r;
                c = R(i,k)/r;
                J = [ c s; -s c];
                R(i:i+1,k:n) = J*R(i:i+1,k:n);  % rotate so that R(i+1, k) is zero
           end
        end
    end
end

