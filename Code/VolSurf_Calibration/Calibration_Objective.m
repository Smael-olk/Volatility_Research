function total_error = Calibration_Objective(params, data)
%CALIBRATION_OBJECTIVE Fast NIG Calibration (No Regularization).
%   INPUTS: params, data
%   Handles Call/Put/Mixed data automatically based on flags.

    total_error = 0;
    
    % Loop through Maturities
    for i = 1:size(data, 1)
        F0 = data{i,3};
        B  = data{i,2};
        Dt = data{i,1};
        MarketPrices = data{i,4};
        Strikes      = data{i,5};
        Flags        = data{i,9}; % 1=Call, 2=Put
        
        if isempty(MarketPrices), continue; end
        
        % --- FAST PRICING (Single FFT) ---
        % 1. Price everything as if it were a CALL.
        % This is efficient because FFT depends on Params+Time, not Strike/Type.
        AllAsCalls = NIG_Pricer(params, F0, Strikes, B, Dt, 1);
        
        % Safety Check
        if any(isnan(AllAsCalls)) || ~isreal(AllAsCalls)
            total_error = 1e10; return;
        end

        ModelPrices = AllAsCalls;
        
        % 2. Convert actual Puts using Parity
        idx_P = (Flags == 2);
        if any(idx_P)
             % Put = Call - B*(F - K)
             ModelPrices(idx_P) = AllAsCalls(idx_P) - B * (F0 - Strikes(idx_P));
             
             % Floor negative prices (numerical noise)
             ModelPrices(idx_P) = max(0, ModelPrices(idx_P));
        end

        % --- ERROR CALCULATION ---
        sq_error = sum((MarketPrices - ModelPrices).^2);
        
        total_error = total_error + sq_error;
    end
end