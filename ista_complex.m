function x = ista_complex(A, b, lambda, max_iter)
% Iterative Shrinkage-Thresholding Algorithm (ISTA) for complex numbers
% A: matrix
% b: vector
% lambda: regularization parameter
% max_iter: maximum number of iterations

% Compute the Lipschitz constant of the gradient
L = norm(A, 2)^2;

% Calculate eta
eta = 1/L;

% Initialize the estimate of the solution vector
x = zeros(size(A, 2), 1);

% Iterate until convergence or maximum number of iterations reached
for iter = 1:max_iter
    % Compute the gradient
    r = b - A*x;
    grad = A'*r; % h
    
    % Update the estimate of the solution vector
    x_new = x + eta*grad;
    x = soft_threshold(x_new, lambda*eta);
    
    % Check for convergence
    if norm(x - x_new) < 1e-6
        break
    end
end
end

function y = soft_threshold(z, T)
% Soft thresholding function for complex numbers
% z: complex number
% T: threshold parameter

y = sign(z) .* max(abs(z) - T, 0);
end