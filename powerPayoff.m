function payoff = powerPayoff(Type, Price, Strike, Power)
% Description of powerPayoff goes here

S = Price;  K = Strike;
m = Power;

switch Type
    case 1
        V = max(S - K, 0).^m;
    case 0
        V = max(K - S, 0).^m;
end

payoff = V;
end