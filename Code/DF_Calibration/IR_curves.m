function [discount_factors,forward_prices,expiries,TradeDate,Spot,mkt] = IR_curves(current_options_file,spot_file)
%IR_CURVES Summary of this function goes here
%   Detailed explanation goes here

    mkt = Option_estraction_hist(current_options_file, spot_file, 0.1, 0.6, 5,50);
    TradeDate = datetime(current_options_file(1:10), 'InputFormat', 'yyyy-MM-dd');
    Spot = mkt.spot; 
    expiries = mkt.datesExpiry;     % final vector for all expiries
    discount_factors = [];     % final vector for all discount factors
    forward_prices = [];    % final vector for all forward prices
    Spot = mkt.spot;
    for i=1:length(expiries)
        Expiry = expiries(i);
        strikes = mkt.strikes(i).value;
        callsmkt = (mkt.callAsk(i).prices + mkt.callBid(i).prices) / 2;
        putsmkt  = (mkt.putAsk(i).prices + mkt.putBid(i).prices) / 2;
        
        syntheticforward = callsmkt - putsmkt;
    
        StrikeLimit = 500;
        keep_mask = (strikes >= (Spot - StrikeLimit)) & (strikes <= (Spot + StrikeLimit));
        strikes = strikes(keep_mask);
        callsmkt = callsmkt(keep_mask);
        putsmkt = putsmkt(keep_mask);
    
        syntheticforward = callsmkt - putsmkt;
        
        Y_data = syntheticforward'; 
        X_data = -strikes';
        X_design = [ones(size(X_data)), X_data];
        B_coeffs = X_design \ Y_data;
        
        discount_factors(i) = B_coeffs(2);
        forward_prices(i) = B_coeffs(1) / B_coeffs(2);
    end
end

