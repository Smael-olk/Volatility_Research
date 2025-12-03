function run_calibration()
%RUN_CALIBRATION Process all files, calibrate NIG model, and save results.

    % --- 1. Setup ---
    data_folder = '.\Dati_Train';
    filelist = dir(fullfile(data_folder, '*.csv')); 
    numFiles = length(filelist);
    spot_file_name = 'spot_SPX_hist.csv'; 
    
    % Storage for results
    Results = struct('Date', {}, 'Params', {}, 'Spot', {}, 'MarketData', {});
    
    % Define the starting guess for the VERY FIRST file
    % (sigma, eta, k)
    
    x0 = [0.2, -0.5, 0.5];
    % --- 3. Optimization Settings ---
    options = optimset('Display', 'iter', 'TolFun', 1e-4, 'MaxFunEvals', 1000);
    

    % --- 2. Loop & Calibrate ---
    for j = 1:numFiles
        % Use fullfile to get the complete path, otherwise IR_curves might not find it
        current_file_path = filelist(j).name;
        
        fprintf('Processing file %d/%d: %s\n', j, numFiles, filelist(j).name);
        
        % Load Data
        [discounts, fwd, expiries, TradeDate, Spot, Option_data] = IR_curves(current_file_path, spot_file_name);
        T_years = years(expiries - TradeDate);
        
        if ~isempty(T_years)
            % Prepare Market Data Cell
            marketdata = cell(length(Option_data.datesExpiry), 8);
            for i = 1:length(Option_data.datesExpiry)
                marketdata{i,1} = T_years(i); % T
                marketdata{i,2} = discounts(i); % B
                marketdata{i,3} = fwd(i); % F0
                marketdata{i,4} = (Option_data.callAsk(i).prices + Option_data.callBid(i).prices)/2; % Price
                marketdata{i,5} = Option_data.strikes(i).value; % Strikes
            end
            
            % --- OPTIMIZATION ---
            
            % --- 4. Run Optimization ---
            disp('Starting Calibration...');
            [optimal_params, fval] = fminsearch(@(x) Calibration_Objective(x, marketdata), x0, options);

            x0= optimal_params;
            
            % --- STORE RESULTS ---
            Results(j).Date = TradeDate;
            Results(j).Spot = Spot;
            Results(j).Params = optimal_params;
            Results(j).T_years = T_years;
            Results(j).Discounts = discounts;
            Results(j).Fwds = fwd;
            
            % --- UPDATE VARIABLES FOR NEXT LOOP ---
            
        end
    end
    
    % --- 3. Save to Disk ---
    save('CalibratedParams.mat', 'Results');
    disp('Calibration complete. Results saved to CalibratedParams.mat');
end