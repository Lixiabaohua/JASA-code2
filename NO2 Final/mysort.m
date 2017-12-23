function SortDat = mysort(mydat )
% resort observation according to time order (day number and hours of day) 
% first sorting 
[sortTime, indTime] = sort(mydat(:,7));  % sort according to 'hours of day'
temp = mydat(:,1:6);                             % response and covariates
sortTemp = temp(indTime,:);                % sort response and covariates 
daynum = mydat(:,8);                           % observation day 
sortDayTime = daynum(indTime);         % sort day number
% second sorting: according to day 
[sortDay, indDay] = sort(sortDayTime);  % sort the sorted day number
SortTimeobse = sortTemp(indDay,:);     % resort response and covariates according to day number
SortTimeDay = sortTime(indDay);         % resort 'hours of day' according to day number

SortDat = [SortTimeobse SortTimeDay sortDay];  %observations in time order
end

