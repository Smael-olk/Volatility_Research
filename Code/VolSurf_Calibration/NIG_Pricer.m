function [Price_model] = NIG_Pricer(params, F0, Strikes, B, Dt, flag)
%NIG_PRICER Prices European options. 
%   Inputs: params, F0, Strikes, B, Dt, flag (1=Call, 2=Put)

    % Unpack parameters
    sigma = params(1);
    eta   = params(2);
    k     = params(3);
    
    % Fixed parameters
    alpha = 1/2; 
    
    % --- CORRECTION HERE ---
    M = 12; % M is the POWER (exponent). N will be 2^12 = 4096 points.
    
    param = 0.001; % Integration step
    
    % --- CONSTRAINT CHECK ---
    if sigma <= 0 || k <= 0
        Price_model = NaN(size(Strikes));
        return;
    end
    
    % --- 1. Calculate Damping Parameter 'a' ---
    p = (1/2 + eta) + sqrt((1/2 + eta)^2 + 2 * (1 - alpha) / (sigma^2 * k));
    a = p / 2; 
    
    % --- 2. Define L(z) and Characteristic Function (phi) ---
    L = @(z, k_val, sig_val) exp(Dt / k_val * (1 - sqrt(1 + 2 * k_val * z * sig_val^2)));
    
    phi = @(u) exp(-1i * u * log(L(eta, k, sigma))) .* ...
               L(((u).^2 + 1i * u * (1 + 2 * eta)) / 2, k, sigma);
               
    % --- 3. Compute Call Price (Base Calculation) ---
    C_model = PriceCallOption(Strikes, F0, B, phi, M, param, a);

    % --- 4. Handle Flag (Call vs Put) ---
    if nargin < 6 || flag == 1
        % Return Call Price
        Price_model = C_model;
        
    elseif flag == 2
        % Return Put Price using Put-Call Parity
        Price_model = C_model - B * (F0 - Strikes);
        
        % Safety constraint for numerical noise
        Price_model(Price_model < 0) = 0; 
    else
        error('Invalid flag. Use 1 for Call, 2 for Put.');
    end
end