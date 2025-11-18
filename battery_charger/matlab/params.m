% default battery_charger params

% Supply
Vrms = 230; 
f = 50;

% Battery
Vbat = 12; 
Rbat = 0.1; 
capacity = 50; 
capUnit = 'Ah';

% SoC & time
SoC_init = 20; 
SoC_target = 80; 
t_charge = [];

% Thyristor
Vt = 0; 
Ileak = 0; 
t_rise = 0; 
t_fall = 0;
