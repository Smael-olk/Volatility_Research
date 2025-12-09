function filtered_mkt = apply_standard_filtering(mkt, discounts, fwd, T_years)
%APPLY_STANDARD_FILTERING Applies the 5 specific option filtering rules.
%   Inputs:
%       mkt: Structure from IR_curves containing raw option data
%       discounts: Discount factors (B)
%       fwd: Forward prices (F)
%       T_years: Time to maturity
%
%   Rules Applied:
%   1. Price < 50 cents REMOVED.
%   2. Spread (Ask-Bid/Mid) > 50% REMOVED.
%   3. Moneyness (K/F) outside 0.90 - 1.10 REMOVED.
%   4. Expiration < 6 days or > 100 days REMOVED.
%   5. No-Arbitrage violations REMOVED.
%   +  Standard OTM Split: Use Puts if K < F, Calls if K >= F.

    % Initialize Output Cell Array (Standard format for calibration)
    % Columns: [Dt, B, F, Price, Strike, ..., ..., ..., Flag]
    filtered_mkt = cell(length(T_years), 9); 
    
    % Loop through each expiry found in the data
    for i = 1:length(T_years)
        Dt = T_years(i);
        B  = discounts(i);
        F0 = fwd(i);
        
        % --- RULE 4: Expiration (6 to 100 days) ---
        days_to_expiry = Dt * 365;
        if days_to_expiry < 6 || days_to_expiry > 100
            continue; % Skip this expiry entirely
        end
        
        % Safety check for index bounds
        if i > length(mkt.strikes), continue; end
        
        % Extract Raw Data Vectors
        strikes  = mkt.strikes(i).value(:);
        call_bid = mkt.callBid(i).prices(:);
        call_ask = mkt.callAsk(i).prices(:);
        put_bid  = mkt.putBid(i).prices(:);
        put_ask  = mkt.putAsk(i).prices(:);
        
        % Calculate Mid Prices and Spreads
        call_mid = (call_bid + call_ask) / 2;
        put_mid  = (put_bid + put_ask) / 2;
        
        call_spread = call_ask - call_bid;
        put_spread  = put_ask - put_bid;
        
        % --- RULE 2: Spread Filter (> 50% of Mid Price) ---
        % (Spread / Mid) <= 0.50
        valid_spread_calls = (call_spread ./ call_mid) <= 0.50;
        valid_spread_puts  = (put_spread ./ put_mid)   <= 0.50;
        
        % --- RULE 3: Moneyness Filter (-10% to +10%) ---
        % 0.90 <= K/F <= 1.10
        moneyness = strikes ./ F0;
        valid_money = (moneyness >= 0.90) & (moneyness <= 1.10);
        
        % --- RULE 5: No-Arbitrage Check ---
        % Call > Max(0, F*B - K*B) -> Call > Intrinsic
        % Put  > Max(0, K*B - F*B) -> Put > Intrinsic
        intrinsic_call = max(0, (F0 - strikes) * B);
        intrinsic_put  = max(0, (strikes - F0) * B);
        
        no_arb_calls = call_mid > intrinsic_call;
        no_arb_puts  = put_mid  > intrinsic_put;
        
        % --- RULE 1: Price Filter (< 50 cents) ---
        price_filter_calls = call_mid >= 0.50;
        price_filter_puts  = put_mid  >= 0.50;
        
        % --- COMBINE: OTM SPLIT + ALL FILTERS ---
        % Logic: 
        % 1. Determine if Strike is PUt or Call side (OTM Split)
        % 2. Apply all boolean masks to that side
        
        % Candidates based on OTM Logic
        is_put_candidate  = (strikes < F0);
        is_call_candidate = (strikes >= F0);
        
        % Final Boolean Masks
        keep_put = is_put_candidate & ...
                   valid_spread_puts & ...
                   valid_money & ...
                   no_arb_puts & ...
                   price_filter_puts;
                   
        keep_call = is_call_candidate & ...
                    valid_spread_calls & ...
                    valid_money & ...
                    no_arb_calls & ...
                    price_filter_calls;
        
        % --- PACK DATA ---
        
        % 1. Puts (Flag = 2)
        idx_p = find(keep_put);
        strikes_p = strikes(idx_p);
        prices_p  = put_mid(idx_p);
        flags_p   = 2 * ones(length(idx_p), 1);
        
        % 2. Calls (Flag = 1)
        idx_c = find(keep_call);
        strikes_c = strikes(idx_c);
        prices_c  = call_mid(idx_c);
        flags_c   = ones(length(idx_c), 1);
        
        % 3. Merge
        final_strikes = [strikes_p; strikes_c];
        final_prices  = [prices_p; prices_c];
        final_flags   = [flags_p; flags_c];
        
        if isempty(final_strikes), continue; end
        
        % 4. Store in Market Data Cell
        filtered_mkt{i, 1} = Dt;
        filtered_mkt{i, 2} = B;
        filtered_mkt{i, 3} = F0;
        filtered_mkt{i, 4} = final_prices;
        filtered_mkt{i, 5} = final_strikes;
        filtered_mkt{i, 9} = final_flags;
    end
end