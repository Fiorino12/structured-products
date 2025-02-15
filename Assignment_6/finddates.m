function [Dates] = finddates(Setdate,increments,flag)
% Finds the payment date given the distance in time
%
% INPUT:
% Setdate:    datenum of today
% increments: increment we want to give in months or years
% flag:       increments in years (flag=0) in months (flag=1)
%
% OUTPUT:
% Dates: business dates of all increments


if (nargin < 3) % number of arguments input
 flag = 0; % default conversion in years
end 

% Initialize today and the vector to store incremented dates
Today = datetime(Setdate, 'ConvertFrom', 'datenum','Format','dd/MM/uuuu');
Dates = datetime(zeros(size(increments)),'ConvertFrom', 'datenum','Format','dd/MM/uuuu');

% Import holidays:
holidays = eurCalendar;
switch flag
    case 0
        % Add to today the increment
        Dates=Today+calyears(increments);
        a = isbusday(Dates, holidays); % check if is a busdate
        % If not find the next bus date 
        Dates(~a)=busdate(Dates(~a),"modifiedfollow",holidays);
        
    case 1
        % Add to today the increment
        Dates=Today+calmonths(increments);
        a = isbusday(Dates, holidays); % check if is a busdate
        % If not find the next bus date 
        Dates(~a)=busdate(Dates(~a),"modifiedfollow",holidays);
        
end

% Convert to datenum
Dates=datenum(Dates);
end