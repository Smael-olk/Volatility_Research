function [Delta, Vega] = Calculate_Greeks_Vanilla(params, S0, K, r, q, T)
%CALCULATE_GREEKS_VANILLA Computes Delta and Vega for a Vanilla Put using NIG.
%   Uses Finite Differences (Bump & Revalue) on the NIG_Pricer.
%
%   INPUTS:
%       params : [sigma, eta, k] (NIG Calibrated Parameters)
%       S0     : Current Spot Price
%       K      : Strike Price of the hedging option
%       r, q   : Risk-free rate and Dividend yield
%       T      : Time to maturity (years)
%
%   OUTPUTS:
%       Delta  : Sensitivity to Spot Price (dP/dS)
%       Vega   : Sensitivity to Volatility (dP/dSigma)

    % --- Configuration ---
    dS_pct = 0.01;      % 1% Spot Perturbation
    dSigma = 0.01;      % 1% Absolute Volatility Perturbation
    
    dS = S0 * dS_pct;
    B = exp(-r * T);    % Discount Factor
    
    % --- 1. Calculate DELTA (Central Difference) ---
    % Note: When Spot (S0) changes, Forward (F0) also changes!
    % F = S * exp((r-q)*T)
    
    % Price Up
    S_up = S0 + dS;
    F_up = S_up * exp((r - q) * T);
    P_up = NIG_Pricer(params, F_up, K, B, T, 2); % Flag 2 = Put
    
    % Price Down
    S_dn = S0 - dS;
    F_dn = S_dn * exp((r - q) * T);
    P_dn = NIG_Pricer(params, F_dn, K, B, T, 2); % Flag 2 = Put
    
    Delta = (P_up - P_dn) / (2 * dS);
    
    % --- 2. Calculate VEGA (Forward Difference) ---
    % Bump Sigma param (params(1))
    
    % Base Price
    F0 = S0 * exp((r - q) * T);
    P0 = NIG_Pricer(params, F0, K, B, T, 2);
    
    % Bumped Price
    params_vega = params;
    params_vega(1) = params(1) + dSigma;
    
    P_vol_up = NIG_Pricer(params_vega, F0, K, B, T, 2);
    
    Vega = (P_vol_up - P0) / dSigma;

end