function [mkt]=Option_estraction_hist(title,namefile_spot,tresh_bid_ask,tresh_penny,DK,tresh_volume)
    %% Call vs Put
    %added threshold for filtering by openinterest
    All_Names=Import_Option_hist(title);
    today = All_Names(2,1);
    All_Names = All_Names(2:end,:);
    
    All_expiry = (table2array(All_Names(2:end,12)));
    Used= ~ismissing(All_expiry);
    
    Names=All_Names(Used,:);
    value_date = table2array(Names(1,1));
    
    expiry=All_expiry(Used);
    
    Put_index=find(table2array(Names(:,end))=='P');
    Call_index=find(table2array(Names(:,end))=='C');
    %discard = [setdiff(table2array(Names(Call_index,7)),table2array(Names(Put_index,7))),setdiff(table2array(Names(Put_index,7)),table2array(Names(Call_index,7)))];
    %Names = Names(~ismember(table2array(Names(:,7)),discard),:);
    %% Data Extraction
    spot_data = import_spot(namefile_spot);
    
    Spot=table2array(spot_data(table2array(spot_data(:,1))==value_date,2));
    
    Strike=Names(:,7);
    
    
    Maturity=unique(expiry);
    All_Call=table2array(Names(Call_index,3:4));
    All_Put=table2array(Names(Put_index,3:4));
    Expiry_call=table2array(Names(Call_index,12));
    Strike_call=table2array(Names(Call_index,7));
    Expiry_put=table2array(Names(Put_index,12));
    Strike_put=table2array(Names(Put_index,7));
    C_vol =table2array(Names(Call_index,6));
    P_vol = table2array(Names(Put_index,6));
    
    dis=0;
    %dist = datenum(Maturity)-datenum(refDate);
    %dist_month = month(Maturity)-month(refDate)+(year(Maturity)-year(refDate))*12;
    %Maturity =Maturity((dist_month<=6| (ismember(month(Maturity),[3,6,9,12])&(year(Maturity)-year(refDate))<=1)|ismember(month(Maturity),[6,12]))&dist>3);
    used=ones(size(Maturity));
    
    for i=1:length(Maturity)
        index_call = (Expiry_call==Maturity(i));
        index_put=(Expiry_put==Maturity(i));
        discard =  [setdiff(Strike_call(index_call),Strike_put(index_put));setdiff(Strike_put(index_put),Strike_call(index_call))];
        Call = All_Call(index_call&(~ismember(Strike_call,discard)),:);
        Put =All_Put(index_put&~ismember(Strike_put,discard),:);
        Strike = Strike_call(index_call&(~ismember(Strike_call,discard)));
    
        index1=find((Call(:,2)~=0).*(Call(:,1)~=0).*(Put(:,1)~=0).*(Put(:,2)~=0));
        index3=index1(find((Call(index1,1)-Call(index1,2))./Call(index1,2)<tresh_bid_ask));
        index4=index3(find((Put(index3,1)-Put(index3,2))./Put(index3,2)<tresh_bid_ask));
        index5=index4(find((C_vol(index4) >= tresh_volume) & (P_vol(index4) >= tresh_volume)));% added filtering by OI
        
        valid_prices=index5(find((C_vol(index5)>tresh_penny*DK).*(P_vol(index5)>tresh_penny*DK)));
        size_s(i) = length(valid_prices);
        if((length(valid_prices)>4)&(datenum(Maturity(i)-table2array(today))>3))
            strikes(i-dis).value=Strike(valid_prices)';
            callAsk(i-dis).prices=Call(valid_prices,1)';
            callBid(i-dis).prices=Call(valid_prices,2)';
            putAsk(i-dis).prices=Put(valid_prices,1)';
            putBid(i-dis).prices=Put(valid_prices,2)';
            Volume_call(i-dis).volume = C_vol(valid_prices);
            Volume_put(i-dis).volume = P_vol(valid_prices);
        
            else
                used(i)=0;
                dis=dis+1;
        end
    end
    mkt.datesExpiry=Maturity(used==1);
    mkt.callBid=callBid;
    mkt.callAsk=callAsk;
    mkt.putAsk=putAsk;
    mkt.putBid=putBid;
    mkt.strikes=strikes;
    mkt.spot=Spot;
    mkt.Volume_call = Volume_call;
    mkt.Volume_put = Volume_put;

end