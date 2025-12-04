function [Delta, Vega] = Calculate_Greeks_NIG(params, S0, r, q, ValuationDate)
%CALCULATE_GREEKS_NIG Computes sensitivities for the Miami Certificate.
%   Uses Finite Differences on the Monte Carlo pricing engine.
%
%   INPUTS:
%       params        : [sigma, eta, k]
%       S0            : Current Spot Price
%       r, q          : Interest rate and Dividend yield
%       ValuationDate : Current Date
%
%   OUTPUTS:
%       Delta : Sensitivity to Spot Price ($ Change per $1 move in Spot)
%       Vega  : Sensitivity to Volatility ($ Change per 100% move in Sigma)

    % --- Configuration ---
    dS_pct = 0.01;      % 1% Spot Perturbation
    dSigma = 0.01;      % 1% Absolute Volatility Perturbation
    
    dS = S0 * dS_pct;
    
    % --- 1. Calculate Baseline Price ---
    [P0, ~] = Price_Miami_Certificate_Function(params, S0, r, q, ValuationDate);
    
    % --- 2. Calculate DELTA (Central Difference) ---
    % Formula: (Price(S+dS) - Price(S-dS)) / (2*dS)
    
    % Price Up
    [P_up, ~] = Price_Miami_Certificate_Function(params, S0 + dS, r, q, ValuationDate);
    
    % Price Down
    [P_down, ~] = Price_Miami_Certificate_Function(params, S0 - dS, r, q, ValuationDate);
    
    Delta = (P_up - P_down) / (2 * dS);
    
    % --- 3. Calculate VEGA (Forward Difference) ---
    % Formula: (Price(Sigma+dSigma) - Price(Sigma)) / dSigma
    % We bump the first parameter of NIG (Sigma)
    
    params_vega = params;
    params_vega(1) = params(1) + dSigma;
    
    [P_vol_up, ~] = Price_Miami_Certificate_Function(params_vega, S0, r, q, ValuationDate);
    
    Vega = (P_vol_up - P0) / dSigma;
    
    % Display for sanity check
    % fprintf('   [Greeks] Delta: %.4f | Vega: %.2f\n', Delta, Vega);

end