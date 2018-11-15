% Numerical Linear Algebra, L.N. Trefethen and D. Bau III
% Lecture 9, Problem 1
%% Lecture code
x = (-128: 128)'/128;  % Set x to a discretization of [- 1,1].
A = [x.^0, x.^1, x.^2, x.^3] ; % Construct Vandermonde matrix.
[Q, R] = qr(A, 0); % Find its reduced QR factorization
scale = Q(257,:);  % Select last row of Q.
Q = Q*diag( 1./scale) ;  %Rescale columns by these numbers.

%% Plotting for (a)
figure
plot(x, Q);
ylim([-1.1, 1.1])
legend('P_0', 'P_1', 'P_2', 'P_3')
xlabel('x')
ylabel('P')
title('Approximate first four Legendre Polynomials')

%% Obtaining and plotting the errors for (b)
P = [x.^0, x, 3/2*x.^2- 1/2, 5/2*x.^3- 3/2*x]; % Set up true values of polynomials
approximationErrors = Q-P;

figure
plot(x, approximationErrors);
legend('P_0', 'P_1', 'P_2', 'P_3')
xlabel('x')
ylabel('Error')
title('Approximation error for first four Legendre Polynomials')

%% Calculate and plot the norms 

error = zeros(3, 4);
error(1,:) = mean(abs(approximationErrors)); % 1-norm
error(2,:) = cellfun(@norm, num2cell(approximationErrors, 1)); %2-norm
error(3,:) = max(abs(approximationErrors)); %infinity-norm
power = [0, 1, 2, 3];
norms = {'1', '2', 'infinity'};

% Plot
figure 
    for i=1:3
       subplot(3, 1, i)
       bar(power,error(i,:));
       suptitle('Approximation error norms')
       title(strcat( norms(i), '-norm')) 
    end
    
    
 %% Track how approximation error changes with the grid becoming finer    
 
 kmin = 8;
 kmax = 25;
 approximationErrorsNorms = zeros(kmax-kmin, 4);
 figure
 hold on
for k=8:25 
    x = [-2^(k): 2^(k)]'/2^(k);  % Set x to a discretization of [- 1,1].
    A = [x.^0, x.^1, x.^2, x.^3] ; % Construct Vandermonde matrix.
    [Q, R] = qr(A, 0); % Find its reduced QR factorization
    scale = Q(2^(k+1)+1,:);  % Select last row of Q.
    Q = Q*diag( 1./scale) ;  %Rescale columns by these numbers.
    P = [x.^0, x, 3/2*x.^2- 1/2, 5/2*x.^3- 3/2*x]; % Set up true values of polynomials
    approximationErrors = Q-P;
    plot(x, approximationErrors);
    approximationErrorsNorms(k-kmin+1,:) = cellfun(@norm, num2cell(approximationErrors, 1)); %2-norm
end
xlabel('x')
ylabel('Error')
title('Approximation errors')

%% Plot 2-norm changing with grid size

figure
plot(kmin:kmax,  approximationErrorsNorms(:,3))
title('Approximation error 2-norm and grid spacing for P_2')
xlabel('Negative of power of 2 in creating the grid')
ylabel('2-norm')

