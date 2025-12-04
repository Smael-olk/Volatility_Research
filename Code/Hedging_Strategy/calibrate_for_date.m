function [optimal_params, S0, r_2y, q_2y, fwd, T_years, Option_data] = calibrate_for_date(date_str)
%CALIBRATE_FOR_DATE Calibrates NIG model for a specific date string.
%   INPUT: date_str (e.g., '2017-12-11')
%   OUTPUT: 
%       optimal_params: [sigma, eta, k]
%       S0: Spot Price
%       r_2y, q_2y: Risk-free rate and Div yield for 2-year maturity (Certificate horizon)
%       fwd, T_years: Full Forward curve and Time vector

    % --- 1. CONFIGURATION ---
    data_folder = '.\Dati_Train';
    spot_file_name = 'spot_SPX_hist.csv'; 
    CalibMode = 'Both'; % Calls + Puts for skew stability
    
    % Certificate Maturity Date (Fixed)
    FinalDate = datetime('09-Dec-2019');
    
    % --- 2. FILE LOADING ---
    % Handle different date formats if necessary, assuming 'yyyy-mm-dd' input
    filelist = dir(fullfile(data_folder, ['*' date_str '*.csv']));
    
    if isempty(filelist)
        error('No option file found for date %s', date_str);
    end
    
    current_file_path = fullfile(filelist(1).folder, filelist(1).name);
    
    % Load Curves
    [discounts, fwd, expiries, TradeDate, Spot, Option_data] = IR_curves(current_file_path, spot_file_name);
    T_years = years(expiries - TradeDate);
    S0 = Spot;
    
    if isempty(T_years)
        error('Data loaded but T_years is empty for %s', date_str);
    end

    % --- 3. DATA FILTERING & PREPARATION ---
    marketdata = cell(length(Option_data.datesExpiry), 9); 
    
    for i = 1:length(Option_data.datesExpiry)
        Dt = T_years(i);
        B = discounts(i);
        F0 = fwd(i);
        
        if Dt <= 0, continue; end
        
        RawStrikes = Option_data.strikes(i).value(:);
        RawCallPrices = (Option_data.callAsk(i).prices(:) + Option_data.callBid(i).prices(:))/2;
        RawPutPrices  = (Option_data.putAsk(i).prices(:)  + Option_data.putBid(i).prices(:))/2;
        
        % Calculate Deltas (Black-76 Proxy)
        sigma_proxy = 0.20;
        d1 = (log(F0 ./ RawStrikes) + (sigma_proxy^2/2)*Dt) / (sigma_proxy*sqrt(Dt));
        Delta_Calls = B * normcdf(d1);
        Delta_Puts  = Delta_Calls - B; 
        
        % Filter: 10% to 90% Delta (Include Wings for Skew)
        idx_calls = find(Delta_Calls >= 0.10 & Delta_Calls <= 0.90);
        idx_puts  = find(abs(Delta_Puts) >= 0.10 & abs(Delta_Puts) <= 0.90);
        
        % Combine based on Mode 'Both'
        MixedStrikes = [RawStrikes(idx_calls); RawStrikes(idx_puts)];
        MixedPrices  = [RawCallPrices(idx_calls); RawPutPrices(idx_puts)];
        TypeFlags    = [ones(length(idx_calls), 1); 2 * ones(length(idx_puts), 1)];
        
        marketdata{i,1}=Dt; marketdata{i,2}=B; marketdata{i,3}=F0;
        marketdata{i,4}=MixedPrices; marketdata{i,5}=MixedStrikes; marketdata{i,9}=TypeFlags; 
    end

    % --- 4. CONSTRAINED OPTIMIZATION (fmincon) ---
    x0 = [0.15, -0.5, 0.5]; 
    
    % Strict Bounds to prevent parameter explosion
    % Sigma > 0, Eta NEGATIVE, K reasonable
    lb = [0.01, -5.00, 0.01]; 
    ub = [1.00, -0.01, 3.00]; 
    
    objective_function = @(x) Calibration_Objective(x, marketdata);
    
    options = optimoptions('fmincon', 'Display', 'off', ... % Silent mode for loop
        'Algorithm', 'sqp', 'MaxFunctionEvaluations', 2000);
        
    [optimal_params, ~] = fmincon(objective_function, x0, [],[],[],[], lb, ub, [], options);
    
    % --- 5. CALCULATE IMPLIED RATES (r, q) FOR CERTIFICATE MATURITY ---
    % Interpolate to find rates specifically for the Certificate's maturity (Dec 2019)
    
    Time_to_Maturity = years(FinalDate - TradeDate);
    
    % 1. Risk Free Rate (r) from Discount Factors
    % r = -ln(B) / T
    % Use linear interpolation on log-discounts
    r_continuous_curve = -log(discounts) ./ T_years;
    r_2y = interp1(T_years, r_continuous_curve, Time_to_Maturity, 'linear', 'extrap');
    
    % 2. Dividend Yield (q) from Forward Prices
    % q = r - ln(F/S) / T
    F_2y = interp1(T_years, fwd, Time_to_Maturity, 'linear', 'extrap');
    q_2y = r_2y - log(F_2y / S0) / Time_to_Maturity;
    
    % Sanity Check for q
    if q_2y > 0.05 || q_2y < -0.01
        warning('Calibrated q=%.2f%% is unrealistic. Capping at 2%%.', q_2y*100);
        q_2y = 0.02;
    end

end