
% Optimal Sizing of Microgrid: Stadium in Norway
%% DENSYS 2022
%% By Emmanuel AF MOMPREMIER



% Cleaning
close all
clear
clc

% Time definition 
% A day will be defined by 24 hours by steps of 30 minutes :

global deltaT Time
deltaT = 0.5;   % half an hour
Time = 0:deltaT:23.5;

% Irradiation power and temperature recovery
A = readmatrix('Data_11_March.xlsx'); %dataset for the day with highest load in winter

Pirr_data = A(1:24,2); %solar irradiance in data (Wh/h)
Temp_data = A(1:24,4); %temperature data (deg C)
Time_data = 0:24; %timescale considered
Pirr_data(25) = Pirr_data(1);
Temp_data(25) = 0.5*(Temp_data(24)+Temp_data(1));
 
% Interpolation and plots
Pirr = interp1(Time_data,Pirr_data,Time);
Temp = interp1(Time_data,Temp_data,Time);
figure(1)
plot(Time,Pirr,'o-','LineWidth',2)
axis([0 24 0 inf],'auto y')
xticks(0:2:24)
xlabel('Day hours')
ylabel('Irradiation Power (W/m^2)')
grid
figure(2)
plot(Time,Temp,'o-','LineWidth',2)
axis([0 24 0 inf],'auto y')
xticks(0:2:24)
xlabel('Day hours')
ylabel('Temperature (°C)')
grid
%% 
% 
% 
% PV power 

global P_PV_STC % PV array rated power in standard test conditions (STC) 
global G_STC    % solar radiation under STC
global CT       % temperature coefficient       The resistance-change factor per degree Celsius of temperature change is called the temperature coefficient of resistance. This factor is represented by the Greek lower-case letter “alpha” (α). A positive coefficient for a material means that its resistance increases with an increase in temperature
global T_STC    % temperature STC
global eta_PV   % panels efficiency (electrical)
global GSR       % global solar radiation
global T_PV     % temperature

P_PV_STC = 380; % W
G_STC = 1000;  % W/m^2
CT = -0.34e-2; % per degree
T_STC = 25;    % celsius degrees
eta_PV = 0.95; 

GSR = Pirr;     % W/m^2
T_PV = Temp;   % celsius degrees
% 
% Plot of one panel produced power during a day

figure(3)
plot(Time,comp_P_PV(1),'LineWidth',2)
axis([0 24 0 inf],'auto y')
xticks(0:2:24)
xlabel('Day hours')
ylabel('One panel power (W)')
grid
% 
% 
% Load Profile recovery and plot

global P_load
P_load_data = 1e3*A(1:24,3);
P_load_data(25) = 0.5*(P_load_data(24)+P_load_data(1))
P_load = interp1(Time_data,P_load_data,Time);

figure(4)
plot(Time,P_load,'o-','LineWidth',2)
axis([0 24 0 inf],'auto y')
xticks(0:2:24)
xlabel('Day hours')
ylabel('Load power profile')
grid
%% 
% 
% 
% Battery and system specifications including investment costs

global C_cell eta_bat eta_DCDC eta_ACDC Cinv_bat Cinv_PV Cinv_cell Cinv_panel C_grid

C_cell = 30*3.6; % cell capacity (W.h) 30 Ah / 3.6 V     
eta_bat = 0.95;  % charging efficiency

eta_DCDC = 0.98; % DC/DC converter efficiency
eta_ACDC = 0.95; % AC DC converter efficiency

Cinv_bat =470 ; % € / kWh  700 470     #we have changed the prices based on data given in a paper studying the energy system in that speciic stadium
Cinv_PV = 7400; % € / kW  %975  7400

Cinv_cell = Cinv_bat*1e-3*C_cell;    %price was for kWh, but cell capacity was for Wh
Cinv_panel = Cinv_PV*1e-3*P_PV_STC;

C_grid = 0.04e-3; % € / kWh

E_load = cumsum(P_load)*deltaT;          % one day energy of the load
E_panel = cumsum(comp_P_PV(1))*deltaT;   % one day energy of one single panel
Npv_min = nearest(E_load/E_panel);       % minimal number of panels
%% 
% 
% 
% Optimal design problem definition
% *Variables :* 
%% 
% * SOC_ini : intial state of battery charge
% * Npv : number of PV panels
% * Ncell : number of battery cells
%% 
% *Objective :*
%% 
% * Cinv : total investissment cost of PV and battery
% * (Optionnal) P_lost : power lost 
%% 
% *Constraints :*
%% 
% * Final SOC > initial SOC
% * SOC min < SOC < SOC max
% * Load always supplied

% X = ( Npv, Ncell, SOC_ini, )

% bounds for SOC constraints
global SOC_min SOC_max
SOC_min = 0.05;
SOC_max = 0.95;
Npv_max = round(5330 / 1.6); %max area covered with panels given in a paper researching that specific stadium:5330 m2 
% area of 1 panel = 1. m^2
Ncell_min = 100;
Ncell_max = 25900;
SOC_ini_min = SOC_min*1e2;
SOC_ini_max = SOC_max*1e2;

% lower and upper bounds constraints

Xl = [Npv_min;Ncell_min;SOC_ini_min];
Xu = [Npv_max;Ncell_max;SOC_ini_max];

obj = @(X) X(1)*Cinv_panel + X(2)*Cinv_cell;  % investment cost

% no linear constraints
% non-linear constraints are defined by a specific function : "constraints"

options = optimoptions(@ga);
[Xopt,cost] = ga(obj,3,[],[],[],[],Xl,Xu,@constraints,[1 2],options)

[P_pv,P_bat,SOC,P_lost,P_unsup] = comp_Power(Xopt);

E_sold = sum(deltaT*P_lost)*365;
rev_min = C_grid*E_sold;

figure(5)
plot(Time,P_load,Time,P_pv,Time,P_bat,Time,P_lost,Time,P_unsup)
legend('Load','PV','Battery','Lost','Unsupplied','Location',"bestoutside")


%% 
% 
%% 
%% Functions definition
% 
% PV power

function P_PV = comp_P_PV(Npv)
    % compute the power produced by Npv panel during a day
    global P_PV_STC G_STC CT T_STC eta_PV GSR T_PV
    P_PV = Npv*eta_PV*P_PV_STC/G_STC*(1+CT*(T_PV-T_STC)).*GSR;
end
%% 
% 
% Powers computation

function [P_pv,P_bat,SOC,P_lost,P_unsup] = comp_Power(X)
    % compute the powers in the different elements
    global deltaT P_load eta_DCDC eta_ACDC eta_bat C_cell SOC_min SOC_max

    Npv = X(1);
    Ncell = X(2);
    SOC_ini = X(3)*1e-2;
    
    P_pv = comp_P_PV(Npv);
    n = length(P_pv);
    P_unsup = zeros(1,n);
    P_bat =  zeros(1,n);
    P_lost = zeros(1,n);
    SOC = zeros(1,n+1);
    SOC(1) = SOC_ini;

    for k = 1:n
        if ( (eta_DCDC*P_pv(k)) >= (P_load(k))/eta_ACDC) 
            if (SOC(k) < SOC_max)
                P_bat(k) = -eta_DCDC*(eta_DCDC*P_pv(k)-P_load(k)/eta_ACDC);
                SOC(k+1) = min(SOC(k)-eta_bat*P_bat(k)*deltaT/(Ncell*C_cell),SOC_max);
            else
                P_lost(k) = eta_DCDC*P_pv(k)-P_load(k)/eta_ACDC; %surplus of energy not stored
                
                SOC(k+1) = SOC(k);
            end
        else
            if (SOC(k) > SOC_min)  
                P_bat(k) = -1/eta_DCDC*(eta_DCDC*P_pv(k)-P_load(k)/eta_ACDC);
                SOC(k+1) = max(SOC(k)-P_bat(k)*deltaT/(Ncell*C_cell),SOC_min);
            else
                P_unsup(k) = -(eta_DCDC*P_pv(k)-P_load(k)/eta_ACDC);
                SOC(k+1) = SOC(k);
            end
        end
    end
end


% Optimization problem constraints
function [c,ceq] = constraints(X)
    % constraints of our optimization problem
    global SOC_min SOC_max
    
    [~,~,SOC,~,P_unsup] = comp_Power(X);
        
    c1 = SOC_min - min(SOC);
    c2 = max(SOC) - SOC_max;
    c3 = max(P_unsup);
    c = [c1 c2 c3]; 
    ceq = [];
end
