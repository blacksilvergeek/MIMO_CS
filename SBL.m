function mu = SBL(Y,A,Delta,convergence,tMax)
% Y = A*mu + n, n~Normal(0,Delta^2)
% convergence is the error threshold, when the difference of adjacent gamma
% is smaller than convergence value, we treat it as convergence.
% tMax is the max number of iteration.
gamma0 = ones(width(A),1);
for i=1:tMax
    % E-step
    variance = (Delta^(-2).*(A'*A)+diag(1./gamma0))^(-1);
    mean = Delta^(-2).*(variance*A'*Y);
    % M-step
    gamma = mean.*mean+diag(variance);
    % check the convergence
    if norm(gamma-gamma0,2)/norm(gamma0,2) < convergence
        break;
    end
    gamma0 = gamma;
end
% compute the final mu.
variance = (Delta^(-2).*(A'*A)+diag(gamma)^(-1))^(-1);
mu = Delta^(-2).*(variance*A'*Y);