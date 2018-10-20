% Numerical Linear Algebra, L.N. Trefethen and D. Bau III
% Lecture 4, Problem 3
function  c4p3(A)
% This function plots the right singular vectors on a unit circle and
% then their image under A
% Input is the 2x2 matrix A
% Function returns no outputs

    % Calculate the right singular vectors
    [V L] = eig(A'*A);
    % Images of v under A
    Us = A*V; 
    % Create the unit circle and its image
    x = [0:0.01:2*pi];
    circ = vertcat(cos(x), sin(x)); % circle
    circA = (A*circ)'; % after applying A to it

    % Plot
    figure

    % First the left singular vectors
    subplot(2,1,1)
    plotv(V) % plot the vectors
    hold on  
    plot(cos(x), sin(x)) % Plot the unit circle
    pbaspect([1 1 1]) % Constrain the aspect ratio so that the circle remains a circle
    % Plot images of u
    subplot(2,1,2)
    hold on
    plotv(Us) % Plot left singular vectors time singular values
    plot(circA(:,1), circA(:,2)) % add the ellipse
    %pbaspect([1 1 1])
    hold off
end