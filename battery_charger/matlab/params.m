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
t_charge = inf;

% Thyristor
alpha_deg = 0:2:180; 
Vt = 10; 
Rth = 0;
Ileak = 10; 
t_rise = 10; 
t_fall = 10;

% Simulation params 
dt = 1/(300*f);

% Visualization
enablePlots = true;
alpha = 30;