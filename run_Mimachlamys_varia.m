close all; 
global pets 

pets = {'Mimachlamys_varia'}; 
check_my_pet(pets); 

estim_options('default'); 
estim_options('max_step_number', 5e2);
estim_options('max_fun_evals', 5e3);

estim_options('pars_init_method', 2); %2 to start from pars_init file, and 1 to start from result mat file
estim_options('results_output', -3); %3 to save everything 
estim_options('method', 'no'); %'nm' pour estimation, 'no' for no estimation

% estim 1.1

estim_pars; 


%%%% 10/06/2023 -- problem of not reaching puberty
% try to let the parameters like this, changed the f value in the predict
% and do not estimate (line 50-51)

% return
% load('results_Mimachlamys_varia.mat');
% other_param = statistics_st('abj', par, C2K(20), 1); %to have other parameters like s_M
% disp('s_M = ')
% disp(other_param.s_M)
% disp('p_Am = ')
% disp(other_param.p_Am)
% disp('E_0 = ')
% disp(other_param.E_0) % current estimation E_0 = 3.10^-3 J cf Eline
% 
% 
return

%% estim 1.2
estim_options('pars_init_method', 1);
% estim_options('results_output', 3);

estim_pars;

% estim_pars;

% estim_pars;
% load('results_Mimachlamys_varia.mat');
% other_param = statistics_st('abj', par, C2K(18), 1); 
% disp('s_M = ')
% disp(other_param.s_M)
% disp('p_Am = ')
% disp(other_param.p_Am)
% disp('E_0 = ')
% disp(other_param.E_0)
% 



%% conversion results.mat to pars_init

%mat2pars_init( )


%%
 %    'results_output':
  %      0     - only saves data to .mat (no printing to html or screen and no figures) - use this for (automatic) continuations 
  %      1, -1 - no saving to .mat file, prints results to html (1) or screen (-1), shows figures but does not save them
  %      2, -2 - saves to .mat file, prints results to html (2) or screen (-2), shows figures but does not save them
  %      3, -3 - like 2 (or -2), but also prints graphs to .png files (default is 3)
  %      4, -4 - like 3 (or -3), but also prints html with implied traits
  %      5, -5 - like 4 (or -4), but includes related species in the implied traits
  %      6     - like 5, but also prints html with population traits