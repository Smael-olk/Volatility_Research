function [S_paths, Times] = NIG_Simulate(params, S0, r, q, T, N_steps, N_sim)
%NIG_SIMULATE_WRAPPER
%   1. Converts your specific [sigma, eta, k] to Canonical [delta, alpha, beta]
%   2. Runs the simulation using the Subordination method.

    % --- 1. Unpack Your Specific Parameters ---
    sigma = params(1); % 0.0861
    eta   = params(2); % 12.4
    k     = params(3); % 0.543

    % --- 2. Perform the Mapping (Derived from your NIG_Pricer) ---
    % Beta (Asymmetry) - Note the negative sign and 0.5 shift
    beta_can = -(eta + 0.5);
    
    % Alpha (Steepness) - Derived from the square root argument match
    alpha_can = sqrt(beta_can^2 + 1 / (k * sigma^2));
    
    % Delta (Scale)
    delta_can = sigma / sqrt(k);



    % --- 3. Run the Subordination Simulation ---
    [S_paths, Times] = NIG_Simulate_Subordination(alpha_can, beta_can, delta_can, S0, r, q, T, N_steps, N_sim);
end

function [S_paths, Times] = NIG_Simulate_Subordination(alpha, beta, delta, S0, r, q, T, N_steps, N_sim)
%   Standard NIG Simulator using Canonical Parameters [alpha, beta, delta]

    dt = T / N_steps;
    Times = linspace(0, T, N_steps + 1);
    
    % Gamma factor (Scaling factor)
    gamma_val = sqrt(alpha^2 - beta^2);
    
    % Inverse Gaussian Parameters for the subordinator
    mu_IG  = delta * dt / gamma_val;
    lam_IG = (delta * dt)^2; 

    % Martingale Drift Correction
    % Ensure alpha > |beta+1| for risk-neutral measure existence
    if alpha <= abs(beta + 1)
        warning('Risk Neutral Constraint Warning: alpha <= |beta+1|');
    end
    
    adj = delta * (gamma_val - sqrt(alpha^2 - (beta + 1)^2));
    drift_rn = (r - q) - adj;
    
    Log_Returns = zeros(N_sim, N_steps);
    
    for i = 1:N_steps
        % A. Generate Random Time Change (Subordinator Y)
        Y = random_inverse_gaussian(mu_IG, lam_IG, N_sim);
        
        % B. Generate Brownian Motion Z
        Z = randn(N_sim, 1);
        
        % C. Construct NIG Increment
        dX = drift_rn * dt + beta * Y + sqrt(Y) .* Z;
        
        Log_Returns(:, i) = dX;
    end
    
    Log_Path = cumsum(Log_Returns, 2);
    S_paths = [ones(N_sim, 1)*S0, S0 * exp(Log_Path)];
end

function x = random_inverse_gaussian(mu, lambda, n)
%   Efficient Sampler for Inverse Gaussian Distribution
    nu = randn(n, 1);
    y = nu.^2;
    mu2 = mu^2;
    term1 = (mu2 .* y) / (2 * lambda);
    term2 = (mu / (2 * lambda)) .* sqrt(4 * mu * lambda .* y + mu2 .* y.^2);
    x = mu + term1 - term2;
    u = rand(n, 1);
    mask = (u > (mu ./ (mu + x)));
    x(mask) = mu2 ./ x(mask);
end