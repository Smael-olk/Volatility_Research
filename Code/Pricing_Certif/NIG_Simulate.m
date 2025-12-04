function [S_paths, Times] = NIG_Simulate(params, S0, r, q, T, N_steps, N_sim)
%NIG_SIMULATE Generates S&P 500 paths using NIG increments.
%   FIXED: Uses small 'a' for CDF calculation and High-Res grid for daily steps.

    dt = T / N_steps;
    Times = linspace(0, T, N_steps + 1);
    
    sigma = params(1);
    eta   = params(2);
    k     = params(3);
    
    % --- 1. Characteristic Function for ONE Step (dt) ---
    L_fun = @(z) exp(dt / k * (1 - sqrt(1 + 2 * k * z * sigma^2)));
    
    phi_NIG = @(u) exp(-1i * u * log(L_fun(eta))) .* ...
                   L_fun((u.^2 + 1i * u * (1 + 2 * eta)) / 2);
    
    % --- 2. Drift Correction ---
    phi_at_minus_i = phi_NIG(-1i);
    martingale_correction = (r - q)*dt - log(phi_at_minus_i);
    
    phi_sim = @(u) exp(1i * u * martingale_correction) .* phi_NIG(u);
    
    % --- 3. Grid Settings for Daily Simulation ---
    % Daily returns are small (~1%). We need a fine grid dx (~0.05%).
    M = 14; 
    N = 2^M;
    target_dx = 0.0005; 
    param = 2 * pi / (N * target_dx); 
    
    % --- 4. Generate Random Increments ---
    a = 1/2; 
    
    Total_Samples = N_sim * N_steps;
    
    % Generate log-returns
    X_all = SimulateFromCF(phi_sim, M, param, a, Total_Samples);
    
    % Reshape into paths
    X_matrix = reshape(X_all, [N_sim, N_steps]);
    
    % --- 5. Build Paths ---
    Log_Path = cumsum(X_matrix, 2);
    S_paths = [ones(N_sim, 1)*S0, S0 * exp(Log_Path)];
    
end