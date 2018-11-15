function Q = formQ(W)
% formQ This function takes output of house() and produces
% the matrix Q in the QR factorization of A based on alg. 10.3
% Inputs: mxn matrix W housing v vectors from Householder algorithm
% Outputs: orthogonal matrix Q

    % Extract size of A and create output matrix
    [m, n] = size(W);
    Q = eye(m); % use e_1...e_m as x vectors in inialization
    % Main loop
    for j=1:m % go over columns
        for k=n:-1:1 % create column
            Q(k:m,j) = Q(k:m,j) - 2*W(k:m,k)*(W(k:m,k)'*Q(k:m,j));
        end
    end
end
