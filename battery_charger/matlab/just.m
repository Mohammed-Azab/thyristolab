if ~exist('var', 'var')
    var = 0; 
end

run params.m

switch var
    case 0
        full_wave_ct_charger(Vrms, f, Vbat, Rbat, capacity, capUnit);
    case 1
        full_wave_ct_charger(Vrms, f, Vbat, Rbat, capacity, capUnit, ...
            't_charge', t_charge, 'SoC_init', SoC_init); 
    case 2
        full_wave_ct_charger(Vrms, f, Vbat, Rbat, capacity, capUnit, ...
            'Vt',Vt ,'Ileak', Ileak, 't_rise', t_rise,'t_fall',t_fall); 
    otherwise
        error('inputNum must be 0, 1, or 2');
end

