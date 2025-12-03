function  dates=business_day(dates,Weekend)
hol=holidays(dates-1000,dates+1000);
if(isbusday(dates,hol,Weekend)~=1)
    	if( lbusdate(year(dates),month(dates),hol,Weekend)<dates)
         dates=lbusdate(year(dates),month(dates),hol,Weekend);
        else
        	dates=busdate(dates,1,hol,Weekend);
        end
     end