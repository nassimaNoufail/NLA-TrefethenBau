% Numerical Linear Algebra, L.N. Trefethen and D. Bau III
% Lecture 9, Problem 3

% Note: matrix A supplied externally
[U, SingV, V] = svd(A);
%% Plot
figure 
subplot(1, 2, 1)
plot(x, singV)
xlim([1 15])
xlabel('Value number')
ylabel('\sigma')
subplot(1, 2, 2)
semilogy(x,singV)
xlim([1 15])
xlabel('Value number')
ylabel('log \sigma')
suptitle('Singular values of matrix A')

%% Construct approximations 

Anu = zeros([size(A), rank(A)]);
Vv = V';
Anu(:, :, 1) = singV(1)*U(:,1)*V(:, 1)';

for i=2:rank(A)
    Anu(:,:,i) = Anu(:,:,i-1) + singV(i)*U(:,i)*V(:, i)';
end

% Plotting loop
figure
for j=1:rank(A)
   subplot(5, 2, j)
   pcolor(Anu(:,:,j))
   colormap(gray)
   set(gca,'Ydir','reverse')
   title(['Rank ' num2str(j) ' approximation'])
end