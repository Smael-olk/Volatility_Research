function [C_model] = NIG_Pricer(params, F0, Strikes, B, Dt)
    % Unpack the parameters that the optimizer will vary
    sigma = params(1);
    eta   = params(2);
    k     = params(3);
    
    % Fixed parameters
    alpha = 1/2; 
    M = 15;
    param = 0.001; % deltaK (Integration step)

    % --- CONSTRAINT CHECK ---
    % NIG parameters must be positive (except eta). .
    if sigma <= 0 || k <= 0
        C_model = NaN(size(Strikes));
        return;
    end

    %  Calculate Damping Parameter 'a' 
    p = (1/2 + eta) + sqrt((1/2 + eta)^2 + 2 * (1 - alpha) / (sigma^2 * k));
    a = p / 2; 

    % Define L(z) (NIG Laplace Exponent Component) 
    L = @(z, k_val, sig_val) exp(Dt / k_val * (1 - sqrt(1 + 2 * k_val * z * sig_val^2)));

    % --- 3. Define Characteristic Function (phi) ---
    phi = @(u) exp(-1i * u * log(L(eta, k, sigma))) .* ...
               L(((u).^2 + 1i * u * (1 + 2 * eta)) / 2, k, sigma);

    % --- 4. Call the Core Lewis Integration Function ---
    C_model = PriceCallOption(Strikes, F0, B, phi, M, param, a);
end