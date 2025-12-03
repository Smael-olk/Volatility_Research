function [discount_factors,forward_prices,expiries,TradeDate,Spot,mkt] = IR_curves(current_options_file,spot_file)
%IR_CURVES Extracts market data and calibrates Discount Factors/Forward prices.

    % 1. Load Data (Requires Full Path)
    % We pass the full path to the extraction function
    mkt = Option_estraction_hist(current_options_file, spot_file, 0.1, 0.6, 5, 50);
    
    % 2. Extract Trade Date from Filename
    % First, strip the path to get just the filename (e.g., 'options_2017-06-08.csv')
    [~, fname, ~] = fileparts(current_options_file);
    
    % Use Regex to find the date pattern (Supports 'yyyy-MM-dd' or 'yyyyMMdd')
    % Looks for "20" followed by digits
    date_str_dash = regexp(fname, '20\d{2}-\d{2}-\d{2}', 'match', 'once'); % 2017-06-08
    date_str_compact = regexp(fname, '20\d{6}', 'match', 'once');          % 20170608
    
    if ~isempty(date_str_dash)
        TradeDate = datetime(date_str_dash, 'InputFormat', 'yyyy-MM-dd');
    elseif ~isempty(date_str_compact)
        TradeDate = datetime(date_str_compact, 'InputFormat', 'yyyyMMdd');
    else
        % Fallback: Try the fixed index if regex fails (assuming old format)
        try
            TradeDate = datetime(fname(1:10), 'InputFormat', 'yyyy-MM-dd');
        catch
            error('Could not extract a valid date from filename: %s', fname);
        end
    end

    % 3. Initialize Variables
    Spot = mkt.spot; 
    expiries = mkt.datesExpiry;    
    discount_factors = zeros(length(expiries), 1); % Pre-allocate for speed
    forward_prices = zeros(length(expiries), 1);
    
    % 4. Loop to calibrate B and F
    for i=1:length(expiries)
        Expiry = expiries(i);
        strikes = mkt.strikes(i).value;
        callsmkt = (mkt.callAsk(i).prices + mkt.callBid(i).prices) / 2;
        putsmkt  = (mkt.putAsk(i).prices + mkt.putBid(i).prices) / 2;
        
        % Put-Call Parity: C - P = B * (F - K)
        % C - P = (B*F) - (B)*K
        % Y     = Intercept + Slope * X
        % Let Y = C - P
        % Let X = -K
        % Slope = B (Discount Factor)
        % Intercept = B*F
        
        % Filter for liquid strikes around spot (StrikeLimit)
        StrikeLimit = 500;
        keep_mask = (strikes >= (Spot - StrikeLimit)) & (strikes <= (Spot + StrikeLimit));
        
        % Apply Filter
        strikes_filt = strikes(keep_mask);
        callsmkt_filt = callsmkt(keep_mask);
        putsmkt_filt = putsmkt(keep_mask);
    
        syntheticforward = callsmkt_filt - putsmkt_filt;
        
        Y_data = syntheticforward'; 
        X_data = -strikes_filt';
        
        % Linear Regression: Y = Intercept + Slope * X
        X_design = [ones(size(X_data)), X_data];
        
        % Robust check: Ensure enough points for regression
        if length(Y_data) < 2
            % Fallback if filter removed too many points
            discount_factors(i) = 1; 
            forward_prices(i) = Spot;
            continue;
        end

        B_coeffs = X_design \ Y_data;
        
        % Extract Results
        discount_factor_i = B_coeffs(2); % Slope is B
        term_i = B_coeffs(1);            % Intercept is B*F
        
        discount_factors(i) = discount_factor_i;
        forward_prices(i) = term_i / discount_factor_i;
    end
end