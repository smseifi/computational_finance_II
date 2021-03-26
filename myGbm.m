function pathMat = myGbm(m, s_zero, mu, sigma, T, n)

timeStep= T/n;

pathMat = zeros(m, n+1);
pathMat (:, 1) = s_zero;
for i = 1:m
    phiVect = randn(n, 1);
    rhs = (mu - sigma^2/2) * timeStep + sigma * phiVect * sqrt(timeStep);
    pathMat(i, 2:end) = exp(cumsum(rhs)+log(s_zero));
end