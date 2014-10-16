%% Create time-series files for SCOPE run
%  Data are from 1. Barn tower data by Andrew Richardson
%                2. Weather data by Harvard Forest LTER
%  SCOPE need the following files:
%                1. a time vector: t_.dat
%                2. a year vector: year_.dat
%                3. incoming shortwave radiation: Rin_.dat
%                4. incoming longwave radiation:  Rli_.dat
%                5. air pressure:  p_.dat
%                6. air temperature above the canopy: Ta_.dat
%                7. Vapor pressure above the canopy: ea_.dat
%                8. wind speed: u_.dat
%  In addition, if you have the priori information on Vcmax or total chl,
%  you can put them into the following two arrays (and store in the file):
%                1. vcmax: table_Vcmax_.dat
%                2. chl:   table_Cab_.dat
%  Author: Xi Yang (geoxiyang@gmail.com)
%  History: v1.0 Oct.13, 2014

clear all
clc

%% 1. Read input file

[data_30min,txt,raw] = xlsread('/Volumes/XiYangResearch/Projects/9.Fluorescence/7.WeatherData/Harvard Barn Radiometric Data Master 12 Feb 2014.xlsx' ...
                        ,'B27534:DB33773'); 

[data2,txt2,raw2]    = xlsread('/Volumes/XiYangResearch/Projects/9.Fluorescence/7.WeatherData/hf001-10-15min-m.xlsx' ...
                        ,'A296737:AC309217');
% [str,remain]         = strtok(txt2(:,1),'T');
% rem                  = strtok(remain,'T');



%% 2. Arrange arrays for output

%  2.1 time
startday = 170;
endday   = 299;
tvec     = startday:(30/(24*60)):(endday+1-30/(24*60));
tvec     = tvec';

%  2.2 year
year     = 2013;
yvec     = zeros(numel(tvec),1);
yvec(:,1)= year;

%  2.3 shortwave radiation
rinvec   = data_30min(:,5);

%  2.4 longwave radiation
rlivec   = data_30min(:,11);

%  2.5 air pressure
pvec     = data2(1:2:12479,16);

%  2.6 air temperature
Tavec    = data_30min(:,104);

%  2.7 vapor pressure above the canopy
%  this is calculated using relative humidity and air temperature
%  according to Monteith & Norman Intro to env. biophy. 2013 2nd Ed.
RHvec    = data_30min(:,105);
eavec    = (6.11*10.^((17.5*Tavec)./(240.97+Tavec))).*(RHvec/100);

%  2.8 wind speed
uvec     = data2(1:2:12479,18);

%% 3. Optional table for vcmax and chl

% %  3.1 Vcmax: first row DOY, second row number
% vcmax(:,1) = [170,190,210,240,270];
% vcmax(:,2) = [90,80,60,40,30];
% 
% %  3.2 chl: first row DOY, second row number
% chl(:,1)   = [170,190,210,240,270];
% chl(:,2)   = [50,40,30,20,10];


%% 4. output files

output_dir = '/Volumes/XiYangResearch/src/SCOPE/data/input/dataset HF_ts/';

dlmwrite([output_dir 't_.dat'],tvec);
dlmwrite([output_dir 'year_.dat'],yvec);
dlmwrite([output_dir 'Rin_.dat'],rinvec);
dlmwrite([output_dir 'Rli_.dat'],rlivec);
dlmwrite([output_dir 'p_.dat'],pvec);
dlmwrite([output_dir 'Ta_.dat'],Tavec);
dlmwrite([output_dir 'ea_.dat'],eavec);
dlmwrite([output_dir 'u_.dat'],uvec);

% optional output
% dlmwrite([output_dir 'table_Vcmax_.dat'],vcmax,'\t');
% dlmwrite([output_dir 'table_Cab_.dat'],chl,'\t');
















