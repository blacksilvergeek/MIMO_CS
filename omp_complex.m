function x = omp_complex(A, b, k)
% Orthogonal Matching Pursuit (OMP) algorithm for complex numbers
% A: matrix
% b: vector
% k: sparsity level

[~, n] = size(A);
x = zeros(n, 1);
residual = b;
support = [];
for i = 1:k
    corr = abs(A'*residual);
    [~, j] = max(corr);
    support = [support; j];
    % reduce support to unique elements
    support = unique(support);
    

    Aj = A(:, support);
    xj = pinv(Aj)*b;
    x(support) = xj;
    residual = b - Aj*xj;
    % if get 4 elements, break
    if length(support) == 4
        break;
    end
end
end































