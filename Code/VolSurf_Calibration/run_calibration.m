function run_calibration()
%RUN_CALIBRATION Process files, calibrate NIG model using Standard Filtering.
%   Features: Uses apply_standard_filtering. 
%   NOW INCLUDES: RMSE and MAPE reporting on Price Errors.

    % --- 1. CONFIGURATION ---
    data_folder = '.\Dati_Train';
    spot_file_name = 'spot_SPX_hist.csv'; 
    
    % --- SETUP ---
    filelist = dir(fullfile(data_folder, '*.csv')); 
    numFiles = length(filelist);
    
    % Storage for results (Added RMSE and MAPE fields)
    Results = struct('Date', {}, 'Params', {}, 'Spot', {}, ...
                     'RMSE', {}, 'MAPE', {}, ... % <--- NEW STORAGE
                     'MarketData', {});
    
    % Initial guess for Day 1 [alpha, beta, delta]
    current_guess = [0.2, -0.5, 0.5]; 
    
    options = optimset('Display', 'off', 'TolFun', 1e-4, 'MaxFunEvals', 1000);
    disp(['Starting calibration for ' num2str(numFiles) ' files using Standard Filtering.']);
    
    % --- 2. LOOP & CALIBRATE ---
    for j = 1:numFiles
        current_file_path = fullfile(filelist(j).folder, filelist(j).name);
        fprintf('Processing file %d/%d: %s\n', j, numFiles, filelist(j).name);
        
        % 1. Extract Raw Data and Curves
        [discounts, fwd, expiries, TradeDate, Spot, Option_data] = IR_curves(current_file_path, spot_file_name, 0.1, 0.6, 5, 10);
        T_years = years(expiries - TradeDate);
        
        if ~isempty(T_years)
            
            % 2. APPLY STANDARD FILTERING
            marketdata = apply_standard_filtering(Option_data, discounts, fwd, T_years);
            
            % Check if valid data remains
            if all(cellfun(@isempty, marketdata(:,1)))
                fprintf('  -> Warning: No valid options left after filtering for this date.\n');
                continue;
            end
            
            % --- OPTIMIZATION ---
            objective_function = @(x) Calibration_Objective(x, marketdata);
            
            [optimal_params, fval] = fminsearch(objective_function, current_guess, options);
            
            % --- NEW: CALCULATE RMSE & MAPE ---
            [rmse_val, mape_val] = calculate_pricing_errors(optimal_params, marketdata);
            
            % --- STORE RESULTS ---
            Results(j).Date = TradeDate;
            Results(j).Spot = Spot;
            Results(j).Params = optimal_params;
            Results(j).T_years = T_years;
            Results(j).Discounts = discounts;
            Results(j).Fwds = fwd;
            
            % Store Error Metrics
            Results(j).RMSE = rmse_val;
            Results(j).MAPE = mape_val;
            
            % Keep Hot Start for next date
            current_guess = optimal_params;
            
            % Updated Print Statement
            fprintf('  -> Done. SSE: %.2f | RMSE: %.4f | MAPE: %.2f%%\n', ...
                    fval, rmse_val, mape_val * 100);
        end
    end
    
    % --- 3. Save to Disk ---
    save('./params/CalibratedParams.mat', 'Results');
    disp('Calibration complete.');
end

% --- HELPER FUNCTION FOR ERROR CALCULATION ---
function [rmse, mape] = calculate_pricing_errors(params, data)
%CALCULATE_PRICING_ERRORS Computes RMSE and MAPE on prices.
    all_market_prices = [];
    all_model_prices  = [];
    
    % Iterate over all maturities (rows in marketdata)
    for i = 1:size(data, 1)
        F0 = data{i,3};
        B  = data{i,2};
        Dt = data{i,1};
        MarketPrices = data{i,4};
        Strikes      = data{i,5};
        Flags        = data{i,9};
        
        if isempty(MarketPrices), continue; end
        
        % 1. Price everything as Calls first (Fast FFT)
        AllAsCalls = NIG_Pricer(params, F0, Strikes, B, Dt, 1);
        
        % 2. Convert to Puts where necessary
        ModelPrices = AllAsCalls;
        idx_P = (Flags == 2);
        if any(idx_P)
             % Put = Call - B*(F - K)
             ModelPrices(idx_P) = AllAsCalls(idx_P) - B * (F0 - Strikes(idx_P));
             ModelPrices(idx_P) = max(0, ModelPrices(idx_P));
        end
        
        % 3. Accumulate arrays
        all_market_prices = [all_market_prices; MarketPrices];
        all_model_prices  = [all_model_prices; ModelPrices];
    end
    
    % --- COMPUTE METRICS ---
    residuals = all_model_prices - all_market_prices;
    
    % RMSE: Root Mean Square Error
    rmse = sqrt(mean(residuals.^2));
    
    % MAPE: Mean Absolute Percentage Error
    % Avoid division by zero by setting a tiny floor for market price
    safe_market_prices = all_market_prices;
    safe_market_prices(safe_market_prices < 1e-6) = 1e-6; 
    
    mape = mean(abs(residuals ./ safe_market_prices));
end