function sigma = Black76Inverse(C_target, F, K, B, T, tol, max_iter)
% Finds Implied Volatility (sigma) via Newton-Raphson.
    
    % Initialize
    sigma = 0.2; % Initial guess (often 0.2 to 0.4 works well)
    tol = 1e-6;
    max_iter = 100;
    
    for i = 1:max_iter
        % 1. Calculate Price and Vega at current sigma
        C_model = Black76Price(F, K, B, T, sigma);
        
        % 2. Calculate Vega (Derivative of Price w.r.t sigma)
        d1 = (log(F / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T));
        N_prime_d1 = (1 / sqrt(2*pi)) * exp(-d1^2 / 2); % Standard Normal PDF
        Vega = B * F * sqrt(T) * N_prime_d1; % Black-76 Vega
        
        % 3. Calculate Error and Check Tolerance
        error = C_model - C_target;
        if abs(error) < tol
            return;
        end
        
        % 4. Newton-Raphson Update: sigma(n+1) = sigma(n) - f(sigma)/f'(sigma)
        sigma = sigma - error / Vega;
        
        % Bound sigma to prevent errors (e.g., sigma > 0)
        sigma = max(1e-6, sigma); 
    end
    
    % If max iterations reached, return NaN or current estimate
    if i == max_iter
        warning('IV solver failed to converge.');
    end
end