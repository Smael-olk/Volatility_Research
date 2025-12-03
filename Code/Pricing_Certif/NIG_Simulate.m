function [S_paths, Times] = NIG_Simulate(params, S0, r, q, T, N_steps, N_sim)
%NIG_SIMULATE Generates S&P 500 paths using NIG increments.
%   Uses the exact CF definition from your calibration.

    dt = T / N_steps;
    Times = linspace(0, T, N_steps + 1);
    
    % Unpack Parameters
    sigma = params(1);
    eta   = params(2);
    k     = params(3);
    
    % --- 1. Define Characteristic Function for ONE Step (dt) ---
    % This matches your NIG_Pricer logic exactly.
    
    % The component L used in your pricer
    % Note: In calibration L depended on Dt. Here it depends on dt (step size).
    L_fun = @(z) exp(dt / k * (1 - sqrt(1 + 2 * k * z * sigma^2)));
    
    % The Characteristic Function of the log-return over dt
    % phi(u) = E[exp(i*u*X)]
    phi_NIG = @(u) exp(-1i * u * log(L_fun(eta))) .* ...
                   L_fun((u.^2 + 1i * u * (1 + 2 * eta)) / 2);
    
    % --- 2. Drift Correction (Risk-Neutral) ---
    % We need E[S(t+dt)] = S(t) * exp((r-q)*dt)
    % Equivalently: E[exp(X)] = exp((r-q)*dt) => phi(-1i) = exp((r-q)*dt)
    % The raw NIG increment might not have this mean. We add a drift 'mu'.
    % phi_corrected(u) = exp(1i*u*mu) * phi_NIG(u)
    % condition: exp(mu + log(phi_NIG(-1i))) = exp((r-q)*dt)
    % mu = (r - q)*dt - log(phi_NIG(-1i))
    
    martingale_correction = (r - q)*dt - log(phi_NIG(-1i));
    
    % Final CF for the simulation
    phi_sim = @(u) exp(1i * u * martingale_correction) .* phi_NIG(u);
    
    % --- 3. Generate Random Increments ---
    % We need N_sim * N_steps random numbers.
    % Optimization: Call SimulateFromCF once for the total batch.
    
    Total_Samples = N_sim * N_steps;
    
    % Settings for FFT Inversion
    M = 12; param = 0.005; a = -0.5;
    
    % Generate ALL log-returns at once (Fastest method)
    X_all = SimulateFromCF(phi_sim, M, param, a, Total_Samples);
    
    % Reshape into paths [N_sim x N_steps]
    X_matrix = reshape(X_all, [N_sim, N_steps]);
    
    % --- 4. Build Paths ---
    Log_Path = cumsum(X_matrix, 2);
    
    % Convert to Price Paths (Add S0 at start)
    S_paths = [ones(N_sim, 1)*S0, S0 * exp(Log_Path)];
    
end