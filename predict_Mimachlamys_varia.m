function [prdData, info] = predict_Mimachlamys_varia(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); %compound parameters = function to computes compound parameters from primary parameters which are defined in par
  vars_pull(par); % unpacks variables from the structure where the primary parameters are
  vars_pull(cPar);  % unpacks variables from the structure where the compound parameters are
  vars_pull(data); % unpacks variables from the structure where the data are
  vars_pull(auxData); % unpacks variables from the structure where the auxiliary data are
  
  % compute temperature correction factors
  % for each data with temperature given
  TC_ab = tempcorr(temp.ab, T_ref, T_A);
  TC_aj = tempcorr(temp.aj, T_ref, T_A);
  TC_tp = tempcorr(temp.tp, T_ref, T_A);  %ajout 24/01/2023
  TC_am = tempcorr(temp.am, T_ref, T_A);
  %TC_Ri = tempcorr(temp.Ri, T_ref, T_A);
  TC_RL = tempcorr(temp.R_L, T_ref, T_A);
  TC_ajnat = tempcorr(temp.aj_nat, T_ref, T_A);


%   TC_19SA = tempcorr(temp.tL19SA, T_ref, T_A); %temperature correction for the data sets tL19SA and LW19SA
%   TC_LW19SA = tempcorr(temp.LW19SA, T_ref, T_A);
  TC_tLlarvae = tempcorr(temp.tLlarvae, T_ref, T_A);
  
  %%% filter for del_M
   if del_M_postmetH >= 1 || del_M_larv >= 1 
    info = 0; prdData = []; return;
   end

%    if f_tL19SA > 3 || f_Tinduff > 1
%    if f_Tinduff < 0.2 || f_tL19SA > 2.5
   if f_tL19SA > 2
    info = 0; prdData = []; return;
   end

  
  %%% 08/02/2023 - update for the del_M with several values depending on the life stage

  % zero-variate data

  % as zero-variate data for Lp Li are issued from L. Regnier-Brisson's in
  % situ monitoring, we will use f_tL19SA for f
  
 % to calculate E_0, we consider in situ f values not the f for the larvae.
 % 
   % initial
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  U_E0 = initial_scaled_reserve(f_tL19SA, pars_UE0); % d.cm^2, initial scaled reserve
  E_0 = J_E_Am * U_E0 * mu_E; % J, initial reserve for kap and v constant
 
   f = f_Tinduff;% Hatchery for early life stages
%   f = f_tL19SA; % try with the other f for maturity problem
%   TC_natural = tempcorr(C2K(14), T_ref, T_A); %temperature in natural environment to see the problem with maturity
    
  % life cycle
  pars_tj = [g k l_T v_Hb v_Hj v_Hp];
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
  
  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Lw_b = L_b/ del_M_larv;                % cm, shell height at birth at f
  aT_b = t_b/ k_M/ TC_ab;           % d, age at birth at f and T
%   aT_b = t_b/ k_M/ TC_natural;    

  % metamorphosis
  L_j = L_m * l_j;                  % cm, structural length at metam
  Lw_j = L_j/ del_M_postmetH;                % cm, shell length at metam at f
  tT_j = t_j / k_M/ TC_aj;             % d, time since fertilisation at metam 
  
%   tT_j = t_j / k_M/ TC_natural;
  
   f = f_tL19SA;% 2019 Sainte Anne
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);

  % metamorphosis
  L_j_nat = L_m * l_j;                  % cm, structural length at metam
  Lw_j_nat = L_j_nat/ del_M_postmetH;                % cm, shell length at metam at f
  tT_jnat = t_j / k_M/ TC_ajnat;             % d, time since fertilisation at metam 
  
  % puberty 
  L_p = L_m * l_p;                  % cm, structural length at puberty at f
  Lw_p = L_p/ del_M_postmetH;                % cm, shell height at puberty at f
  tT_p = (t_p - t_b)/ k_M/ TC_tp;   % d, time since birth at puberty at f and T
  
  % ultimate
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_M_postmetH  ;              % cm, ultimate total shell height at f
  Wd_i = d_V * L_i.^3  *(1 + f * w );  %in g, dry weight without gonad at f

 
  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/TC_am;                % d, mean life span at T
  
   % reproduction (from P. maximus)
  f = f_repro;% reproduction conditions (mix between Tinduff hatchery and Saint-Anne monitoring)
  
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector at T
  RT_L = TC_RL * reprod_rate_j(L_R, f, pars_R);                 % #/d, ultimate reproduction rate at T

  E_R = RT_L * 365/2 * E_0 / kap_R;
  W_at_LR = d_V * (L_R * del_M_postmetH).^3  *(1 + f * w );  %in g, dry weight without gonad at f
  W_ER_at_LR = w_E / mu_E * E_R ;
  GSI = W_ER_at_LR / (W_at_LR + W_ER_at_LR);
  
  % taken from Pinna nobilis but needs checking
  %GSI = 365 * TC_GSI * k_M * g/ F^3/ (F + kap * g * y_V_E);
  %GSI = GSI * ((1 - kap) * F^3 - k_M^2 * g^2 * k_J * U_Hp/ v^2/ s_M^3);
  
  

  % pack to output
  prdData.ab = aT_b;    % age at birth
  prdData.aj = tT_j;    % age at metamorphosis since fertilisation
  prdData.aj_nat = tT_jnat;
  prdData.tp = tT_p;    % d, time since birth at puberty at f and T     
  prdData.am = aT_m;    % life span
  prdData.Lb = Lw_b;    % length at birth
  prdData.Lj = Lw_j;    % length at metamorphosis
  prdData.Lj_nat = Lw_j_nat;
  prdData.Lp = Lw_p;    % cm, shell height at puberty at f              
  prdData.Li = Lw_i;    % ultimate length
  prdData.Wdi = Wd_i;    % ultimate dry weight
%   prdData.Ri = RT_i;    % ultimate reproduction rate
%   
  prdData.R_L = RT_L;    % reproduction rate at L_R
  prdData.GSI = GSI;
  prdData.E_0 = E_0;
  
  % uni-variate data
  
  %%% Time-Length of larvae and juveniles, after birth and before puberty
  %------------------------------------------------------------------------
  % Data from Tinduff hatchery
  % Estimation routine from Engraulis encrasicolus predict
  % modified on the 09/02/2023 to separate larvae from juveniles in the
  % data set 
  
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_Tinduff); %f_Tinduff, the functional response for all data sets coming from Tinduff hatchery
  
  %%% Larvae 
  rT_j = TC_tLlarvae * rho_j * k_M;         % exponential growth rate, corrected with the temperature
  time_bj = tLlarvae(:,1);                  % days, time between birth and metamorphosis
  Lw_b = l_b * L_m/ del_M_larv;             % cm, physical length at birth
  
  ELw_bj = Lw_b * exp(time_bj * rT_j/3);    % estimated length between birth and metamorphosis
  
  %%% Juveniles
  rT_B = TC_tLlarvae * rho_B * k_M;         % von Bertalanffy growth rate, corrected with the temperature
  time_ji = tLjuve(:,1);                    % days, time after metamorphosis
  aT_j = t_j/ k_M/ TC_tLlarvae;             % time at metamorphosis, since fertilisation
  Lw_j = l_j * L_m/ del_M_postmetH;          % cm, physical height at metamorphosis %06/06/23 
  Lw_i = l_i * L_m/ del_M_postmetH;          % cm, physical ultimate height
  ELw_ji = Lw_i - (Lw_i - Lw_j) * exp( - rT_B * (time_ji - aT_j));  % estimated height after metamorphosis
  
  
  %------------------------------------------------------------------------
 
  %%% Time-Length-Ash free dry weight, of juveniles and adults
  %------------------------------------------------------------------------
  % Saint-Anne/Roscanvel data set, from Laure Regnier-Brisson
  % Temperature values from environmental data set, defined here with
  % Fourier transform function and parameters.
  
  % L_metam is the same as the L_j zero-variate data, but defined in the parameters set

  init_cond_tLW = [L_metam*par.del_M_postmetH, ((par.z*par.p_M/par.kap)/par.v)*(L_metam*par.del_M_postmetH)^3, par.E_Hj, 0]; %[L_0, E_0, EH_0, ER_0] %try to start for metamorphosis time
  %del_M_postmetH corresponds to shape coeff for individuals after metamorphosis with shell height measurements
  % tL19SA(1,2) the initial length 
  
 [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_tL19SA);
  L_b = l_b * L_m; 
  L_j = l_j * L_m;
%   a_tL19SA = [a_metam ; tL19SA(:,1)];  %start for metamorphosis time -- a_metam is the same as the a_j zero-variate data, but defined in the parameters set
  
%  [t_sort_19SA, it , it_sort_19SA] = unique(a_tL19SA,'sorted');
 [t_sort_19SA, it , it_sort_19SA] = unique(tL19SA(:,1),'sorted'); %take the time of the data
 
  pars_temp19SA = temp.tL19SA;      % Fourier parameters for temperature (from the environment) : used for tL and tW
  [~, LEHR] = ode45(@ode_LEHR, t_sort_19SA, init_cond_tLW, [], par, f_tL19SA, pars_temp19SA, L_b, L_j); %pars_temp19SA, the Fourier parameters %TC_LW19SA for one temperature value
  L      = LEHR(:,1);               % cm, structural length
  E      = LEHR(:,2);             % J, reserve
%   E_H    = LEHR(:,3);             % J, maturity level
  E_R    = LEHR(:,4);               % J, reproduction buffer
  
  Lw   = L/ del_M_postmetH;          % cm, standard length
  Lw   = Lw(it_sort_19SA);          % reconstruct L
%   ELw_19SA = Lw(2:end);             % mm, standard length
   ELw_19SA = Lw;             % mm, standard length %% pour calcul de L avec valeur initiale
  
  E_R = E_R(it_sort_19SA);          % reconstruct E_R 
%   E_R = E_R(2:end);                 % J, delete the first value which is for the value close to 0

  E = E(it_sort_19SA);              % reconstruct E
%   E = E(2:end);                    % J, delete the first value which is for the value close to 0
   
%   ELw_forW = Lw(2:end);
  EtWw_19SA = (ELw_19SA.*del_M_postmetH).^3 + E_R .* w_E/ mu_E/ d_E + E.* w_E/ mu_E/ d_E;   % g ww, wet weight
  %in this equation we suppose that d_V = d_E; and d_Vd = 1 so it is not
  %written. If we change d_V, for the estimation, we also need to change d_E 
  EtWd_19SA = d_V * EtWw_19SA;                                      % g dw, dry weight

  %length-weight estimation /!\ E_R not considered
  EWd_19SA = d_V * ((LW19SA(:,1)*del_M_postmetH ).^3  *(1 + f_tL19SA * w )) ; %in g, dry weight

  
  %------------------------------------------------------------------------
  
  % pack to output
  prdData.tLlarvae = ELw_bj;
  prdData.tLjuve = ELw_ji;
  prdData.tL19SA = ELw_19SA;  
  prdData.tW19SA = EtWd_19SA;
  prdData.LW19SA = EWd_19SA;

end
  
%%% subfunctions
%--------------------------------------------------------------------------

      function dLEHR = ode_LEHR(t, LEHR, p, f, par_temp, L_b, L_j)
    % function from Engraulis ringens predict: (modified to add E_R)
    % Input: 
    % p: structure 'par' 
    % c: structure 'Cpar' obtained by cPar = parscomp_st(par)
    % f: scaled, scaled functional response, 
    % s_M: scalar, -, acceleration factor post metamorphosis
    % TC, scalar, -, temperature correction factor
    % L_b, scaler, cm, structural length at birth at f
    % L_j, scaler, cm, structural length at metamorphosis at f

    % --------------- unpack LEHR ------------------------------------------
    L   =  max(0,LEHR(1)); % cm, volumetric structural length
    E   =  max(0,LEHR(2)); % J,   energy in reserve 
    EH  =  max(0,LEHR(3)); % J, E_H maturity
    ER  =  max(0,LEHR(4)); % J, E_R reproduction buffer

    % shape correction function:
    if EH < p.E_Hb                              %before birth
        s_M = 1;
    elseif EH >= p.E_Hb && EH < p.E_Hj          %between birth and metamorphosis
        s_M = L/L_b;
    else                                        %after metamorphosis
        s_M = L_j/L_b;
    end
    % Temperature and shape correct the relevant paramters
    T = fnfourier(t, par_temp);     % take the temperature from the Fourier parameters (that were calculated in another script)
    TK = C2K(T);                    % temperature in Kelvin degree
    TC = tempcorr(TK, p.T_ref, p.T_A);  % temperature correction factor
%     TC = par_temp; %test with only one value of temperature correction factor
    vT    = s_M * p.v * TC; 
    pT_Am = s_M * p.z * p.p_M/ p.kap * TC;
    pT_M  = p.p_M * TC; 
    kT_J  = p.k_J * TC; 
    %
    pA  = f * pT_Am * L^2 * (EH >= p.E_Hb);           % J/d, assimilation
    r   = (E * vT/ L - pT_M * L^3/ p.kap)/ (E + p.E_G * L^3/ p.kap);
    pC  = E * (vT/ L - r); % J/d, mobilisation 
    dE  = pA - pC;               % J/d, change in energy in reserve
    dL  = r/ 3 * L;              % cm/d, change in structural length
    dEH = ((1 - p.kap) * pC - kT_J * EH) * (EH < p.E_Hp);    % J/d, change in cum energy invested in maturation (it is implied here that no rejuvenation occurs).
    dER = ((1 - p.kap) * pC - kT_J * p.E_Hp)*(EH >= p.E_Hp); % J/d, change in cum energy invested in reproduction (the maturity maintenance is constant)

    % pack dLEHR
    dLEHR = [dL; dE; dEH; dER];    
    end
    %----------------------------------------------------------------------