%% Q3. part (c)

clc
clear variables 

% ----------------------------------------------------------------------- % 

Type = 1;
S_zero = 1e2;    K = 1e2;
r = 2e-2;    sigma = 23e-2;
T = 1;
m = 1;

N = 250;

[~, sMat, deltMat] = binomialDeltaPowerCall(Type, S_zero, K, r, sigma, ...
    T, m, N);

S = linspace(80, 140, 100);

binomDeltVect = deltMat(:, 1+N/2);
binomPrices = sMat(:, 1+N/2);

interpDeltVect = interpDelta(binomDeltVect, binomPrices, S);

blsDeltVect = zeros(1, 100);
for i = 1:100
    blsDeltVect(1, i) = blsdelta(S(1, i), K, r, 0.5, sigma);
end

f1 = figure(1);
[default_position, ~] = plset;
f1.InnerPosition = default_position;
f1.OuterPosition = default_position;

plot(S, blsDeltVect, 'k-', 'LineWidth', 3.5, 'DisplayName', ...
    'Black--Scholes Delta')
hold on 
plot(S, interpDeltVect, 'r:', 'LineWidth', 2, 'DisplayName', ...
    'Binomial Interpolated Delta')

ax = gca;   ax.FontSize = 15;   ax.LineWidth = 1;

ax.XGrid = 'on';    ax.YGrid = 'on';
ax.XMinorTick = 'on';   ax.XMinorGrid = 'on';
ax.YMinorTick = 'on';   ax.YMinorGrid = 'on';

legend('Location', 'southeast', 'Interpreter', 'LaTex', 'LineWidth', ...
    1.25, 'FontSize', 14);

%% Q3. part (d) + (e)

clc
clear variables 

% ----------------------------------------------------------------------- %

Type = 1;
S_zero = 1e+2;  K = 1e+2;
r = 2e-2;   sigma = 23e-2;  mu = 18e-2;
T = 1;
m = 1;

N = 250;    timeStep = T/N;
M = 62500;

[V_zero, sMat, deltMat] = binomialDeltaPowerCall(Type, S_zero, K, r, ...
    sigma, T, m, N);

pathMat = myGbm(M, S_zero, mu, sigma, T, N);

compound = exp(r * (timeStep));

nVect = 20:20:N;  nl = length(nVect);
interpDeltMat = zeros(M, nl);
interpDeltMat(:, 1) = deltMat(end, 1);
for j = 2:nl
    n = nVect(1, j-1);
    interpDeltMat(:, j) = interpDelta(deltMat(:, n+1), sMat(:, n+1), ...
        pathMat(:, n+1).');
end

deltDiffMat = -diff(interpDeltMat, 1, 2);

bMat = zeros(M, nl);
bMat(:, 1) = V_zero - interpDeltMat(:, 1) .* pathMat(:, 1);

correctionMat = deltDiffMat .* pathMat(:, nVect(1:nl-1)+1);

for j = 2:nl
    bMat(:, j) = bMat(:, j-1) * compound ^ (nVect(1, 1)) + ...
        correctionMat(:, j-1);
end

S_NVect = pathMat(:, end);
V_NVect = powerPayoff(Type, S_NVect, K, m);

Pi_NVect = -V_NVect + interpDeltMat(:, end) .* S_NVect + ...
    compound ^ (nVect(1, 1)) * bMat(:, end);
pNlvect = exp(-r) * Pi_NVect / V_zero;

[VaR, CVaR] = dVaRCVaR(pNlvect, 95e-2);

f2 = figure(2);
[default_position, ~] = plset;
f2.InnerPosition = default_position;
f2.OuterPosition = default_position;

h = histogram(pNlvect, 50, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'k';

ax = gca;   ax.FontSize = 15;   ax.LineWidth = 1;

ax.XGrid = 'off';    ax.YGrid = 'off';
ax.XMinorTick = 'off';   ax.XMinorGrid = 'off';
ax.YMinorTick = 'off';   ax.YMinorGrid = 'off';

disp([VaR, CVaR])

%% Q4. part (d)

clc
clear variables 

% ----------------------------------------------------------------------- %

S_zero = 95;  K = 105;
r = 0;   sigma = 23e-2;  mu = 18e-2;
T = 1;
m = 1;

N = 8e+2;   M = 8e+4;   timeStep = T/N;

compound = exp(r * (timeStep));

nVect = 1:1:N;  nl = length(nVect);

pathMat = myGbm(M, S_zero, mu, sigma, T, N);
C_zeroVect = ones(M, 1) * blsprice(S_zero, K, r, 1, 23e-2);
C_TVect = powerPayoff(1, pathMat(:, end), K, m);

deltMat = zeros(M, nl);
deltMat(:, 1) = 0;
deltMat(:, 2:nl) = max(sign(pathMat(:, nVect(1:nl-1)+1) - K), 0);

deltDiffMat = -diff(deltMat, 1, 2);
correctionMat = deltDiffMat .* pathMat(:, nVect(1, 1:nl-1)+1);

bMat = zeros(M, nl);
bMat(:, 1) = C_zeroVect;
for j = 2:nl
    bMat(:, j) = bMat(:, j-1) * compound ^ (nVect(1, 1)) + ...
        correctionMat(:, j-1);
end

pNlvect = (-C_TVect + deltMat(:, end) .* pathMat(:, end) + ...
    bMat(:, end)) ./ C_zeroVect;

[VaR, CVaR] = dVaRCVaR(pNlvect, 95e-2);

f3 = figure(3);
[default_position, ~] = plset;
f3.InnerPosition = default_position;
f3.OuterPosition = default_position;

h = histogram(pNlvect, 50);
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'k';

ax = gca;   ax.FontSize = 15;   ax.LineWidth = 1;

ax.XGrid = 'off';    ax.YGrid = 'off';
ax.XMinorTick = 'off';   ax.XMinorGrid = 'off';
ax.YMinorTick = 'off';   ax.YMinorGrid = 'off';
disp([mean(pNlvect), var(pNlvect) VaR, CVaR])

%% Q5. part (c)

clc
clear variables 

% ----------------------------------------------------------------------- %

Type = 1;
S_zero = 95;  K = 95;
r = 5e-2;   sigma = 2e-1;
T = 1;
m = 1;

uMu = 32e-2;    dMu = 3e-1;
uP = 4e-1;  dP = 1 - uP;
kappa = uP/(1 - uMu) + dP/(1 + dMu) - 1; 

lambda = 1e-1;

drift = r - sigma^2/2 - lambda * kappa;


N = 8e+2;   M = 25e+3;   timeStep = T/N;
jP = lambda * timeStep;

xVect = log(S_zero) * ones(M, 1);


for j = 1:N
    jMask = rand(M, 1) <= jP;
    isUp = rand(M, 1) <= uP;
    jAmp = isUp .* exprnd(uMu, [M, 1]) + (1 - isUp) .* ...
        -exprnd(dMu, [M, 1]);
    jVect = jMask .* jAmp;
    xVect = xVect + drift * timeStep + sigma * sqrt(timeStep) * ...
        randn(M, 1) + jVect;
end

S_NVect = exp(xVect);
V_NVect = powerPayoff(1, S_NVect, K, m);
V_zero = exp(-r) * mean(V_NVect);

%% Q5. part (d)

kl = 20;
KVect = linspace(70, 120, kl);

impVolVect = zeros(1, kl);
for j = 1:kl
    K = KVect(1, j);
    V_NVect = powerPayoff(1, S_NVect, K, m);
    V_zero = exp(-r) * mean(V_NVect);
    impVolVect(1, j) = blsimpv(S_zero, K, r, 1, V_zero);
end

%%
f4 = figure(4);
[default_position, ~] = plset;
f4.InnerPosition = default_position;
f4.OuterPosition = default_position;

plot(KVect, impVolVect, 'r:', 'LineWidth', 2.5, 'DisplayName', ...
    'Implied Volatility');
hold on
yline(sigma, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Realized Volatility')
ax = gca;   ax.FontSize = 15;   ax.LineWidth = 1;

legend('Location', 'northeast')

ax.XGrid = 'on';    ax.YGrid = 'on';
ax.XMinorTick = 'on';   ax.XMinorGrid = 'on';
ax.YMinorTick = 'on';   ax.YMinorGrid = 'on';

ax.YLim = [.18, .26];

%% Q6. part (a) + (b)

clc
clear variables 

% ----------------------------------------------------------------------- %

load RawData.mat

trainMat = zeros(size(CSTrain, 2), 8);
trainMat(:, 1) = 1;
trainMat(:, 2) = STrain;
trainMat(:, 3) = TauTrain;
trainMat(:, 4) = DeltaTrain;
trainMat(:, 5) = DeltaTrain .* 2;
trainMat(:, 6) = VegaTrain;
trainMat(:, 7) = VegaTrain .* 2;
trainMat(:, 8) = VegaTrain .* DeltaTrain;

testMat = zeros(size(CSTest, 2), 8);
testMat(:, 1) = 1;
testMat(:, 2) = STest;
testMat(:, 3) = TauTest;
testMat(:, 4) = DeltaTest;
testMat(:, 5) = DeltaTest .* 2;
testMat(:, 6) = VegaTest;
testMat(:, 7) = VegaTest .* 2;
testMat(:, 8) = VegaTest .* DeltaTest;

% ----------------------------------------------------------------------- %

cVect = (trainMat) \ (CVTrain ./ CSTrain).';

trainDeltVect = trainMat * cVect;

trainErrorVect = CVTrain.' - (trainDeltVect .* CSTrain.');
trainBench = CVTrain.' - (DeltaTrain .* CSTrain.');

[trainVar, trainCvar] = dVaRCVaR(trainErrorVect, 95e-2);
disp([mean(trainErrorVect), var(trainErrorVect), trainVar, trainCvar])


[trainBenchVar, trainBenchCvar] = dVaRCVaR(trainBench, 95e-2);
disp([mean(trainBench), var(trainBench), trainBenchVar, trainBenchCvar])

% ----------------------------------------------------------------------- %

testDeltVect = testMat * cVect;

testErrorVect = CVTest.' - (testDeltVect .* CSTest.');
testBench = CVTest.' - (DeltaTest .* CSTest.');

[testVar, testCvar] = dVaRCVaR(testErrorVect, 95e-2);
disp([mean(testErrorVect), var(testErrorVect), testVar, testCvar])


[testBenchVar, testBenchCvar] = dVaRCVaR(testBench, 95e-2);
disp([mean(testBench), var(testBench), testBenchVar, testBenchCvar])

f5 = figure(5);
[default_position, ~] = plset;
f5.InnerPosition = default_position;
f5.OuterPosition = default_position;

h = histogram(testErrorVect, 50, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'k';

ax = gca;   ax.FontSize = 15;   ax.LineWidth = 1;

ax.XGrid = 'off';    ax.YGrid = 'off';
ax.XMinorTick = 'off';   ax.XMinorGrid = 'off';
ax.YMinorTick = 'off';   ax.YMinorGrid = 'off';

%% Q6. part (c)

clc
clear variables 

% ----------------------------------------------------------------------- %

load RawData.mat

trainMat = zeros(size(CSTrain, 2), 3);
trainMat(:, 1) = 1;
trainMat(:, 2) = DeltaTrain; 
trainMat(:, 3) = DeltaTrain .^ 2;

testMat = zeros(size(CSTest, 2), 3);
testMat(:, 1) = 1;
testMat(:, 2) = DeltaTest; 
testMat(:, 3) = DeltaTest .^ 2; 

% ----------------------------------------------------------------------- %

trainRHS = ((CVTrain ./ CSTrain).' - DeltaTrain) .* STrain ...
    .* sqrt(TauTrain) ./ (VegaTrain);

cVect = (trainMat) \ trainRHS;

trainDeltVect = DeltaTrain + VegaTrain ./ (STrain .* sqrt(TauTrain))...
    .* (trainMat * cVect);

trainErrorVect = CVTrain.' - (trainDeltVect .* CSTrain.');
trainBench = CVTrain.' - (DeltaTrain .* CSTrain.');

beta = 95e-2; 

[trainVar, trainCvar] = dVaRCVaR(trainErrorVect, beta);
trainV = var(trainErrorVect);
trainMean = mean(trainErrorVect);
disp([trainVar, trainCvar, trainV, trainMean])



[trainBenchVar, trainBenchCvar] = dVaRCVaR(trainBench, beta);
trainBenchV = var(trainBench);
trainBenchMean = mean(trainBench);
disp([trainBenchVar, trainBenchCvar, trainBenchV, trainBenchMean])

% ----------------------------------------------------------------------- %

testDeltVect = DeltaTest + VegaTest ./ (STest .* sqrt(TauTest))...
    .* (testMat * cVect);

testErrorVect = CVTest.' - (testDeltVect .* CSTest.');
testBench = CVTest.' - (DeltaTest .* CSTest.');

[testVar, testCvar] = dVaRCVaR(testErrorVect, beta);
testV = var(testErrorVect);
testMean = mean(testErrorVect);
disp([testVar, testCvar, testV, testMean])

[testBenchVar, testBenchCvar] = dVaRCVaR(testBench, beta);
testBenchV = var(testBench);
testBenchMean = mean(testBench);
disp([testBenchVar, testBenchCvar, testBenchV, testBenchMean])

f6 = figure(6);
[default_position, ~] = plset;
f6.InnerPosition = default_position;
f6.OuterPosition = default_position;

h = histogram(testErrorVect, 50, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'k';

ax = gca;   ax.FontSize = 15;   ax.LineWidth = 1;

ax.XGrid = 'off';    ax.YGrid = 'off';
ax.XMinorTick = 'off';   ax.XMinorGrid = 'off';
ax.YMinorTick = 'off';   ax.YMinorGrid = 'off';





