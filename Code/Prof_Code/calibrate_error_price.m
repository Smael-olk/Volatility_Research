function error_value = calibrate_error_price(params, MarketData, F0, Discounts, T_years)
%CALIBRATE_ERROR Summary of this function goes here
%   Detailed explanation goes here
    sigma_nig = params(1); % Using a different name to avoid conflict
    eta = params(2);
    k = params(3);
    
    alpha = 1/2; 
    total_squared_error = 0;
    
    % --- Loop through all Maturities (T_j) ---
    for j = 1:length(T_years)
        Dt = T_years(j);
        B = Discounts(j);
        
        % Recalculate 'a' for the new parameters
        p = (1/2 + eta) + sqrt((1/2 + eta)^2 + 2 * (1 - alpha) / (sigma_nig^2 * k));
        a = p / 2;
        
        % Define the NIG Characteristic Function for this maturity
        % (Use sigma_nig)
        phi = @(u) exp(-1i * u * log(exp(Dt / k * (1 - sqrt(1 + 2 * k * (u.^2 + 1i * u * (1 + 2 * eta)) / 2 * sigma_nig^2)))));

        % Extract Strikes and Market Prices for this maturity
        Strikes = MarketData.strikes;
        MarketPrices = (MarketData.callAsk(j).value-MarketData.call(j).value) / 2; % <<< Use the observed prices
        
        % 1. Price ALL options using the Lewis Formula (Your existing code)
        % Assuming your function handles the vector of strikes
        C_model_vector = PriceCallOption(Strikes, F0, B, phi, 15, 0.001, a); 

        % 2. Calculate the Squared Price Error
        squared_errors = (MarketPrices - C_model_vector).^2;
        total_squared_error = total_squared_error + sum(squared_errors);
    end
    
    error_value = total_squared_error;
end

