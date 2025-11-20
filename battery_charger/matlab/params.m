% default battery_charger params

% Supply
Vrms = 230; 
f = 50;

% Battery
Vbat = 12; 
Rbat = 0.1; 
capacity = 50; 
capUnit = 'Ah';
SoC_target = 100;

% SoC & time
SoC_init = 20; 
SoC_target = 80; 
t_charge = inf;

% Thyristor
alpha_deg = 0:5:175; 
Vt = 0; 
Rth = 0;
Ileak = 0; 
t_rise = 0; 
t_fall = 0;

% Simulation params 
dt = 1/(300*f);

% Visualization
enablePlots = true;
alpha = 30;