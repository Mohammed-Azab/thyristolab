% default battery_charger params

% Supply
Vrms = 230; 
f = 50;

% Battery
Vbat = 20; 
Rbat = 0.1; 
capacity = 50; 
capUnit = 'Ah';

% SoC & time
SoC_init = 12; 
SoC_target = 100;
t_charging_hours =  0.5; % Simulated User Input
t_charge = inf;

% Thyristor
alpha = 30;
alpha_deg = 0:2:180;         
Vt = 1.5;       
Rth = 0.001;  
Ileak = 0.01;  
t_rise = 1e-6;  
t_fall = 2e-6;  

% Simulation params 
dt = 1/(300*f);

% Visualization
enablePlots = true;
savePlots = false; 
