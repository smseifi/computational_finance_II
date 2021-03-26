function [initialValue, binomialPrices, binomialDeltas] = ...
    binomialDeltaPowerCall(Type, Price, Strike, Rate, Volatility, ...
    Maturity, Power, Periods)
% Description of binomialDeltaPoweCall goes here

s_zero = Price; k = Strike;
r = Rate;
sigma = Volatility;
m = Power;
T = Maturity;   n = Periods;

Delta = T/n;

discount = exp(-r * Delta);

u = exp(sigma * sqrt(Delta) + (r - (sigma^2)/2) * Delta);
d = exp(-sigma * sqrt(Delta) + (r - (sigma^2)/2) * Delta);
q = (exp(r * Delta) - d)/(u - d);

sMat = zeros(n+1, n+1);
for j = 1:n+1
    sMat(n-j+2:n+1, j) = (s_zero * d.^(0:j-1) .* u.^(j-1:-1:0)).';
end

vMat = zeros(n+1, n+1); deltMat = zeros(n, n);
switch Type
    case 1
        vMat(:, n+1) = powerPayoff(1, sMat(:, n+1), k, m);
    case 0
        vMat(:, n+1) = powerPayoff(0, sMat(:, n+1), k, m);
end
for j = n:-1:1
    vMat(n-j+2:n+1, j) = discount * ((1-q)*vMat(n-j+2:n+1, j+1) + ...
        q*vMat(n-j+1:n, j+1));
    
    deltMat(n-j+1:n, j) = (vMat(n-j+1:n, j+1)-vMat(n-j+2:n+1, j+1)) ...
        ./ (sMat(n-j+1:n, j+1)-sMat(n-j+2:n+1, j+1));
end

initialValue = vMat(n+1, 1);
binomialPrices = sMat;
binomialDeltas = deltMat;
end