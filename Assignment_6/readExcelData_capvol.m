function [Data_vol] = readExcelData_capvol( filename, Settlement)
% Reads data from excel
%  It reads the flat volatility surface and relevant dates
%  All input volatilities are in bps units
%
% INPUTS:
%  filename: excel file name where data are stored
%  Settlement: Settlement date in datenum
% 
% OUTPUTS:
%  Data_vol: struct with Settlement, strikes , expiry dates and matrix of volatilities

% Put Settlement in the output
Data_vol.Settlement = Settlement;

% Upload strikes
strikes = xlsread(filename, 1, 'F1:R1');
Data_vol.strikes = strikes/100;

% Upload expiry dates
[~, incr_dates] = xlsread(filename,1,'B2: B17');
incr_dates = str2double(strrep(incr_dates, 'y', ''));
incr_dates = incr_dates([1, 3:end]);
Data_vol.expyear = incr_dates; 
Data_vol.expiries = finddates(Settlement, incr_dates);

% Upload volatilities
volatilities = xlsread(filename,1,'F2:R17');
volatilities = volatilities([1,3:end],:);
Data_vol.flat_volatilities = volatilities*1e-4;

