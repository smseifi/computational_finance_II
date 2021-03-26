function interpDeltas = interpDelta(binomialDeltas, binomialPrices, Prices)
% Description of interpDelta goes here

dVect = [min(binomialDeltas); binomialDeltas]; 
sVect = binomialPrices; 

dVect = dVect(sVect ~= 0);
sVect = sVect(sVect ~= 0);

deltVect = interp1(sVect, dVect, Prices, 'linear', 'extrap');
for j = 1:length(deltVect)
deltVect(1, j) = min(max(deltVect(1, j), min(dVect)), max(dVect));
end

interpDeltas = deltVect;
end

