function run_calibration()
%RUN_CALIBRATION Process files, calibrate NIG model. 
%   Features: ATM Delta Filter, Configurable Mode (Call/Put/Both), NO Regularization.

    % --- 1. CONFIGURATION ---
    data_folder = '.\Dati_Train';
    spot_file_name = 'spot_SPX_hist.csv'; 
    
    % Select Calibration Mode: 'Call', 'Put', or 'Both'
    CalibMode = 'Put'; 
    
    % --- SETUP ---
    filelist = dir(fullfile(data_folder, '*.csv')); 
    numFiles = length(filelist);
    
    % Storage for results
    Results = struct('Date', {}, 'Params', {}, 'Spot', {}, 'MarketData', {});
    
    % Initial guess for Day 1
    current_guess = [0.2, -0.5, 0.5]; 
    
    options = optimset('Display', 'off', 'TolFun', 1e-4, 'MaxFunEvals', 1000);
    disp(['Starting calibration for ' num2str(numFiles) ' files in mode: ' CalibMode]);
    
    % --- 2. LOOP & CALIBRATE ---
    for j = 1:numFiles
        current_file_path = fullfile(filelist(j).folder, filelist(j).name);
        fprintf('Processing file %d/%d: %s\n', j, numFiles, filelist(j).name);
        
        [discounts, fwd, expiries, TradeDate, Spot, Option_data] = IR_curves(current_file_path, spot_file_name);
        T_years = years(expiries - TradeDate);
        
        if ~isempty(T_years)
            marketdata = cell(length(Option_data.datesExpiry), 9); 
            
            % --- DATA PREPARATION ---
            for i = 1:length(Option_data.datesExpiry)
                
                Dt = T_years(i);
                B = discounts(i);
                F0 = fwd(i);
                
                if Dt <= 0, continue; end
                
                % 1. Extract Raw Data (Force Column Vectors)
                RawStrikes = Option_data.strikes(i).value(:);
                RawCallPrices = (Option_data.callAsk(i).prices(:) + Option_data.callBid(i).prices(:))/2;
                RawPutPrices  = (Option_data.putAsk(i).prices(:)  + Option_data.putBid(i).prices(:))/2;
                
                % 2. Calculate Deltas (Black-76) for ATM filtering
                sigma_proxy = 0.20;
                d1 = (log(F0 ./ RawStrikes) + (sigma_proxy^2/2)*Dt) / (sigma_proxy*sqrt(Dt));
                Delta_Calls = B * normcdf(d1);
                Delta_Puts  = Delta_Calls - B; 
                
                % 3. Determine Indices based on ATM Filter
                idx_calls = find(Delta_Calls >= 0.40 & Delta_Calls <= 0.60);
                idx_puts  = find(abs(Delta_Puts) >= 0.40 & abs(Delta_Puts) <= 0.60);
                
                % 4. Select Data Based on CalibMode
                switch CalibMode
                    case 'Call'
                        MixedStrikes = RawStrikes(idx_calls);
                        MixedPrices  = RawCallPrices(idx_calls);
                        TypeFlags    = ones(length(idx_calls), 1); % 1=Call
                        
                    case 'Put'
                        MixedStrikes = RawStrikes(idx_puts);
                        MixedPrices  = RawPutPrices(idx_puts);
                        TypeFlags    = 2 * ones(length(idx_puts), 1); % 2=Put
                        
                    case 'Both'
                        MixedStrikes = [RawStrikes(idx_calls); RawStrikes(idx_puts)];
                        MixedPrices  = [RawCallPrices(idx_calls); RawPutPrices(idx_puts)];
                        TypeFlags    = [ones(length(idx_calls), 1); 2 * ones(length(idx_puts), 1)];
                end
                
                % 5. Store in marketdata
                marketdata{i,1} = Dt; 
                marketdata{i,2} = B; 
                marketdata{i,3} = F0; 
                marketdata{i,4} = MixedPrices; 
                marketdata{i,5} = MixedStrikes;
                marketdata{i,9} = TypeFlags; 
            end
            
            % --- OPTIMIZATION (No Regularization Inputs) ---
            objective_function = @(x) Calibration_Objective(x, marketdata);
            
            [optimal_params, fval] = fminsearch(objective_function, current_guess, options);
            
            % --- STORE RESULTS ---
            Results(j).Date = TradeDate;
            Results(j).Spot = Spot;
            Results(j).Params = optimal_params;
            Results(j).T_years = T_years;
            Results(j).Discounts = discounts;
            Results(j).Fwds = fwd;
            
            % Keep Hot Start
            current_guess = optimal_params;
        end
    end
    
    % --- 3. Save to Disk ---
    save('CalibratedParams.mat', 'Results');
    disp('Calibration complete.');
end