function C = Black76Price(F, K, B, T, sigma)
% Calculates the European Call Price using the Black-76 model.
    d1 = (log(F / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T));
    d2 = d1 - sigma * sqrt(T);
    
    N_d1 = 0.5 * (1 + erf(d1 / sqrt(2))); % Standard Normal CDF
    N_d2 = 0.5 * (1 + erf(d2 / sqrt(2)));
    
    C = B * (F * N_d1 - K * N_d2);
end