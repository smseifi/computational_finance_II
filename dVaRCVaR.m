function [var, cvar] = dVaRCVaR(pNlVect, beta)

sortedpNlVect = sort(pNlVect);

M = length(pNlVect);
idx = floor((1-beta) * M);

var = sortedpNlVect(idx);
cvar = (1/idx) * sum(sortedpNlVect(1:idx));
end