function [par, metaPar, txtPar] = pars_init_Mimachlamys_varia(metaData)

metaPar.model = 'abj'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.T_A = 10000;      free.T_A   = 0;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 
par.z = 0.60618;      free.z     = 1;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 6.5;        free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.027654;     free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.92542;    free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 38.3388;    free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 2343.5335;  free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 1.082e-04; free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hj = 7.116e-04; free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam'; 
par.E_Hp = 3.930e+01; free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 1.565e-08;  free.h_a   = 0;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.0001;     free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.L_R = 5.5;        free.L_R   = 0;   units.L_R = 'cm';         label.L_R = 'length at reproduction  data'; 
par.L_metam = 0.0215;  free.L_metam = 0;   units.L_metam = 'cm';     label.L_metam = 'shell height at metamorphosis for starting data'; 
par.a_metam = 18;     free.a_metam = 0;   units.a_metam = 'd';      label.a_metam = 'age since fertilisation at metamorphosis for starting data'; 
par.del_M_larv = 0.79114;  free.del_M_larv = 0;   units.del_M_larv = '-';   label.del_M_larv = 'shape coefficient for larvae data'; 
par.del_M_postmetH = 0.38339;  free.del_M_postmetH = 0;   units.del_M_postmetH = '-';  label.del_M_postmetH = 'shape coefficient for individuals after metamorphosis, shell height'; 
par.f = 1;            free.f     = 0;   units.f = '-';            label.f = 'scaled functional response for 0-var data'; 
par.f_Tinduff = 0.14155;  free.f_Tinduff = 1;   units.f_Tinduff = '-';    label.f_Tinduff = 'scaled functional response for tLlarvae and tLjuve data from Tinduff Hatchery'; 
par.f_repro = 0.2753;  free.f_repro = 1;   units.f_repro = '-';      label.f_repro = 'scaled functional response for reproduction data'; 
par.f_tL19SA = 1.9997;  free.f_tL19SA = 1;   units.f_tL19SA = '-';     label.f_tL19SA = 'scaled functional response for tL19SA data'; 

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
