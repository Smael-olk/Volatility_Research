function total_error = Calibration_Objective(params, data)
    % params = [sigma, eta, k]
    
    total_error = 0;
    penalty_weight = 1000; % penalty for arbitrage
    
    % Loop through Maturities (if data has multiple)
    for i = 1:length(data{1})
        F0 = data{i,3};
        B  = data{i,2};
        Dt = data{i,1};
        Strikes = data{i,5};
        MarketPrices = data{i,4}; % market peices
        
        ModelPrices = NIG_Pricer(params, F0, Strikes, B, Dt);
        
        % Nan handle
        if any(isnan(ModelPrices))
            total_error = 1e10; 
            return;
        end
        
        % MSE
        sq_error = sum((MarketPrices - ModelPrices).^2);
        
        % 3. ARBITRAGE PENALTY (Butterfly Check)
        if length(ModelPrices) > 2
            Butterfly = ModelPrices(1:end-2) - 2*ModelPrices(2:end-1) + ModelPrices(3:end);
            
            if any(Butterfly < -1e-6) % If negative probability detected
                sq_error = sq_error + penalty_weight;
            end
        end
        
        total_error = total_error + sq_error;
    end
end