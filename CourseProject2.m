%Benjamin Henry
%CEE586: Physical Hydrology
%Course Project 2
%
%
%
%
%This script reads in all the necessary data provided in the assignment
%description and performs all the necessary calculations needed to come up
%with a water balance for the East River Watershed. Additionally, this
%script also plots the SWE in the watershed, which corresponds to task 6 in
%the Course Project deliverables, and creates all other necessary plots.
%
%An Important note is that since all of the data is either given in terms
%of daily or hourly values, the overall water balance produced here will be
%on a daily frequency. Any data given hourly will be averaged on a daily
%scale. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%WATER BALANCE ASSUMPTION

%Assuming that the water balance can be expressed by the following equation
%
% dS(subsurface)/dt + dS(surface)/dt = P - ET - R - G 
%
%In this project we are neglecting groundwater (G = 0), and assuming that
%the discharge rates are a good proxy for the runoff rates (R = discharge)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Reading in all of the data

%SNOTEL SITES (Contains Precipitation and Water Table Accumulation
%(Amount of Water in Snowpack))
SNOTEL1 = readtable('380_STAND_WATERYEAR=2018.csv');
SNOTEL2 = readtable('737_STAND_WATERYEAR=2018.csv');

%FLUX TOWER (Contains Daily Evaporation-Transpiration Values)
FLUX = readtable('flux_tower_et.csv');

%PUMP HOUSE DISCHARGE (Contains Discharge Rates in m^3/s in a 10 min
%resolution)
PUMPHOUSE = readtable('PumpHouse_discharge_wy2018.csv');

%WATER LEVEL (Contains The height level of the underground water table,
%which can be used to infer subsurface storage)
WATERLEVEL = readtable('waterlevel_PLM1_PLM6_soil_wy2018.csv',...
    'NumHeaderLines',14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PREPARING THE DATA (Mainly conversion of units to m^3/day)

%CREATING A TIME ARRAY
Time = table2array(SNOTEL1(:,2));
Time(366) = []; %extra entry in array needs to be removed

%AREA OF WATERSHED
Area = 85; %km^2
Area = Area*1000000; %km^2 to m^2

%PRECIPITATION
%Find total precipitation in inches (total precip is
%taken by adding together the precip measured at both SNOTEL sites)

Tot_Precip = table2array(SNOTEL1(:,5)) + table2array(SNOTEL2(:,5)); 
Tot_Precip = Tot_Precip*0.0254*Area; %inches to meters^3 conversion

%Converting Total Precipitaion into Precipitation Rates by taking the 
%differences between each day's measurements.
for i=1:length(Tot_Precip)-1
    Precip(i,1) = Tot_Precip(i+1)-Tot_Precip(i); %Precipitation Rate in m^3/day
end

%EVAPORATION/TRANSPIRATION
Evap = table2array(FLUX(:,2)); %ET in mm/day
Evap = Evap*0.001*Area; %mm/day to m^3/day

%There were a lot of days where the flux was not measured. In order to
%account for this, a mean ET was taken using nanmean(), which computes the
%mean of all the values while ignoring the NaN values. This mean value was
%used to replace the NaN values.
Mean_Evap = nanmean(Evap);
for i=1:length(Evap)
    if isnan(Evap(i))
        Evap(i) = Mean_Evap; %Getting rid of NaN values
    end
end


%RUNOFF
R = table2array(PUMPHOUSE(:,2));  %Runoff in m^3/s
%Find Runoff Rates in Units of m^3/day. The Runoff Rates are given
%instantaneously in 10 minute intervals or 1 measurement per day, so in 
%order to find the correct units, these these rates must be multiplied by 
%the correct time duration to represent Runoff per Day Values. The 
%time/frequency of measurements changes on day to day which is 
%frustrating to calculate but I try my best to account for this and 
% still produce daily averages 

j = 1; %index for each day
Runoff = zeros(365,1); %creating the array to store runoff values

for i=1:length(R) %for the total amount of measurements
    %this time interval has measurements every 10 minutes but is missing
    %two measurements
    if i < 143
        Runoff(j) = Runoff(j) + R(i)*10*60;
    end
    %this time interval has measurements every 10 minutes
    if i > 143 && i < 9359
        if rem(i,144) == 0
            j = j+1;
        end
        Runoff(j) = Runoff(j) + R(i)*10*60; 
    end
    %this time interval has only 1 measurement per day
    if i > 9359 && i < 9485
        j = j+1;
        Runoff(j) = R(i)*86400;
    end
    %this time interval has measurements every 10 minutes
    if i > 9485
        if rem(i,144) == 0
            j = j+1;
        end
        Runoff(j) = Runoff(j) + R(i)*10*60; 
    end
end

%The above for loop creates daily Runoff/Discharge Values for 364 days. To
%create a runoff rate for the very last day, I simply took an average of
%the daily discharge rates.
Runoff(365) = sum(Runoff)/364;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FIRST ESTIMATE OF WATER BALANCE USING PRECIPITATION, EVAPO-TRANSPIRATION,
%AND RUNOFF

%Net Balance = P - E- R

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WATERBALANCE1_Volume = Precip - Runoff - Evap;
WATERBALANCE1_Meters = WATERBALANCE1_Volume/Area; 

%PLOTTING

%Volumetric Water Balance
figure(1)
plot(Time, WATERBALANCE1_Volume)
title('East River Watershed Water Balance(Water Balance 1)')
ylabel('Change in Storage (m^3/day)')
xlabel('Date')

%Water Balance
figure(2)
plot(Time, WATERBALANCE1_Meters,'Color','r')
title('East River Watershed Water Balance (m/day)(Water Balance 1)')
ylabel('Change in Storage (m/day)')
xlabel('Date')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SECOND ESTIMATE OF WATER BALANCE USING CHANGES IN SUBSURFACE WATER LEVEL 
%ELEVATION AND WATER TABLE HEIGHT EQUIVALENT OF SNOW

%The main assumption is that changes in these storage terms are equivalent
%to the balance using precipitation, evaporation, and runoff. It also 
%assumes that any changes in these storage terms gets factored into 
%evaporation or runoff. This second approximation is used to compare to 
%the first estimation. 

%Net Balance = dS(surface)/dt + dS(subsurface)/dt 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ABOVE SURFACE STORAGE(SNOWPACK)
%EQUIVALENT TO THE PRECIPITATION CALCULATIONS

%Find total water equivalent in snowpack in inches (total water eq. is
%taken by adding together the both SNOTEL sites measurements)
TOT_WTEQ_Snow = table2array(SNOTEL1(:,4)) + table2array(SNOTEL2(:,4)); 
TOT_WTEQ_Snow = TOT_WTEQ_Snow*0.0254*Area; %inches to meters^3 conversion
TOT_WTEQ_Snow(366) = []; %Data is repeated/corrupted and was removed 

%Converting Total Water Eq. into Change in Above Surface Rates by taking the 
%differences between each day's measurements.
for i=1:length(TOT_WTEQ_Snow)-1
    %Change in storage in m^3/day
    WTEQ_Snow(i,1) = TOT_WTEQ_Snow(i+1)-TOT_WTEQ_Snow(i);  
end
WTEQ_Snow(365) = sum(WTEQ_Snow)/364;

%SUBSURFACE STORAGE(WATER TABLE LEVEL)
TOT_WATERLEVEL = table2array(WATERLEVEL(:,2)) + table2array(WATERLEVEL(:,5));
TOT_WATERLEVEL = TOT_WATERLEVEL*Area*1; %m/hr to m^3

%Converting from total total water level to change in subsurface storage
%rates by taking the differences beetween each hours's measurements.
for i=1:length(TOT_WATERLEVEL)-1
    %Change in storage in m^3/hour
    W(i,1) = TOT_WATERLEVEL(i+1)-TOT_WATERLEVEL(i); 
end

l = 1; %index for days
WL = zeros(365,1); %creating the array to store Water Level values

%Convert the hourly rate into a daily rate
for i = 1:length(TOT_WATERLEVEL)-1
    WL(l) = WL(l)+W(i);
    if rem(i,24) == 0
        l = l+1;
    end
end

WATERBALANCE2_Volume = WL + WTEQ_Snow;
WATERBALANCE2_Meters = WATERBALANCE2_Volume/Area;

%PLOTTING

%Volumetric Water Balance
figure(3)
plot(Time, WATERBALANCE2_Volume)
title('East River Watershed Water Balance (Water Balance 2)')
ylabel('Change in Storage (m^3/day)')
xlabel('Date')

%Water Balance
figure(4)
plot(Time, WATERBALANCE2_Meters,'Color','r')
title('East River Watershed Water Balance (m/day)(Water Balance 2)')
ylabel('Change in Storage (m/day)')
xlabel('Date')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%THIRD ESTIMATE OF WATER BALANCE BY COMBINING ALL TERMS

% 0 = P - E - R - dS(surface) - dS(subsurface)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WATERBALANCE3_Volume = WATERBALANCE1_Volume - WATERBALANCE2_Volume;
WATERBALANCE3_Meters = WATERBALANCE3_Volume/Area;

%Volumetric Water Balance
figure(5)
plot(Time, WATERBALANCE3_Volume)
title('East River Watershed Water Balance (Water Balance 3)')
ylabel('Change in Storage (m^3/day)')
xlabel('Date')

%Water Balance
figure(6)
plot(Time, WATERBALANCE3_Meters,'Color','r')
title('East River Watershed Water Balance (m/day)(Water Balance 3)')
ylabel('Change in Storage (m/day)')
xlabel('Date')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculating SWE and Plotting

%This is accommplished by taking the water equivalent data from both SNOTEL
%sites (WTEQ values) and combining them

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculation
SWE = table2array(SNOTEL1(:,4))+ table2array(SNOTEL2(:,4)); %SWE (inches)
SWE = SWE*2.54; %inches to centimeters
SWE(366) = []; %Removing corrupted data entry at the end of the dataset

%Plotting 
figure(7)
plot(Time, SWE)
title('Snow Water Equivalent (SWE)')
ylabel('Height (cm)')
xlabel('Date')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MORE PLOTTING FOR ANALYSIS IN WRITEUP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Evap
figure(8)
subplot(2,1,1)
plot(Time, Evap)
title('Evapotranspiration')
ylabel('Flux (m^3/day)')
xlabel('Date')

%Runoff
subplot(2,1,2)
plot(Time, Runoff)
title('Runoff')
ylabel('Flux (m^3/day)')
xlabel('Date')

%Precip 
figure(9)
plot(Time, Precip)
title('Precipitation')
ylabel('Flux (m^3/day)')
xlabel('Date')

%SubSurfaceStorage Water Level
figure(10)
plot(Time, WL/Area)
title('Change in Water Level (Subsurface Storage)')
ylabel('Flux (m/day)')
xlabel('Date')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


