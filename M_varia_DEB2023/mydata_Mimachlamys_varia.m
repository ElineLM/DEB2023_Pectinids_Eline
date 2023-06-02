function [data, auxData, metaData, txtData, weights] = mydata_Mimachlamys_varia
%% set metaData
metaData.phylum     = 'Mollusca'; 
metaData.class      = 'Bivalvia'; 
metaData.order      = 'Ostreoida'; 
metaData.family     = 'Pectinidae';
metaData.species    = 'Mimachlamys_varia'; 
metaData.species_en = 'Variegated scallop'; 

metaData.ecoCode.climate = {'MC'};
metaData.ecoCode.ecozone = {'MANE','MAE','MAm'};
metaData.ecoCode.habitat = {'0jMp', 'jiMb'};
metaData.ecoCode.embryo  = {'Mp'};
metaData.ecoCode.migrate = {};
metaData.ecoCode.food    = {'biPp'};
metaData.ecoCode.gender  = {'Hsb'}; %bidirectional hermaphrodite
metaData.ecoCode.reprod  = {'O'};

metaData.T_typical  = C2K(16); % K. body temp
metaData.data_0     = {'ab'; 'aj'; 'tp'; 'am'; 'Lb'; 'Lj'; 'Lp'; 'Li';'Ri'};      
metaData.data_1     = {'t-L'; 'L-Wd'}; 

metaData.COMPLETE = 2.5; % using criteria of LikaKear2011       %donnees complètes ou pas suiviant nombre et type donnees

metaData.author   = {'Le Moan & RegnierBrisson'};    
metaData.date_subm = [2023 01 11];              
metaData.email    = {'bas.kooijman@vu.nl'};            
metaData.address  = {'VU University Amsterdam'};   

metaData.curator     = {'Starrlight Augustine'};
metaData.email_cur   = {'sta@akvaplan.niva.no'}; 
metaData.date_acc    = [2021 02 21]; 

%% set data
% zero-variate data

data.ab = 2;      units.ab = 'd';    label.ab = 'age at birth';             bibkey.ab = 'TinduffHatchery';   
  temp.ab = C2K(19.5);  units.temp.ab = 'K'; label.temp.ab = 'temperature';
  comment.ab = 'Based on Tinduff hatchery protocol, with controlled temperature';
data.aj = 18;    units.aj = 'd';    label.aj = 'time since fertilisation at metam'; bibkey.aj = 'TinduffHatchery';
  temp.aj = C2K(18);  units.temp.aj = 'K'; label.temp.aj = 'temperature';
  comment.aj = 'Based on Tinduff hatchery protocol, with controlled temperature';
data.tp = 270;     units.tp = 'd';    label.tp = 'time since birth at puberty'; bibkey.tp = {'RegnierBrisson'};  %ajout 24/01/2023 ==> puberté
  temp.tp = C2K(14);  units.temp.tp = 'K'; label.temp.tp = 'temperature'; 
data.am = 7*365;   units.am = 'd';    label.am = 'life span';                bibkey.am = 'sealifebase';   
  temp.am = C2K(16);  units.temp.am = 'K'; label.temp.am = 'temperature'; 
  comment.am = 'Agree with the original data set of AmP, from Sealifebase website';


  
data.Lb  = 0.0105;   units.Lb  = 'cm';  label.Lb  = 'shell length at birth';   bibkey.Lb  = 'TinduffHatchery';
  comment.Lb = 'Based on Tinduff hatchery protocol, with controlled temperature';
  temp.Lb = C2K(19.5); units.temp.Lb = 'K'; label.temp.Lb = 'temperature';
data.Lj  = 0.0215;   units.Lj  = 'cm';  label.Lj  = 'metamorphosis shell length';   bibkey.Lj  = 'TinduffHatchery';
  comment.Lj = 'Based on Tinduff hatchery protocol, with controlled temperature';
  temp.Lj = C2K(18); units.temp.Lj = 'K'; label.temp.Lj = 'temperature';
data.Lp  = 1.5;    units.Lp  = 'cm';  label.Lp  = 'shell height at puberty'; bibkey.Lp  = 'RegnierBrisson';
  comment.Lp = 'Suivi N20';  %ajout 24/01/2023
  
  %added on 15th May 2023
data.Li = 7.2; units.Li = 'cm'; label.Li = 'infinite shell height'; bibkey.Li = 'RegnierBrisson';
    comment.Li = 'guess from in situ monitoring and Perodou & Latrouite 1981';

data.R_L  = 2*4000000/365; units.R_L  = '#/d'; label.R_L  = 'reprod rate at a large size';     bibkey.R_L  = 'TinduffHatchery';   %test *2 for the Ri, 27/01/2023 -- impact of the kappa, but not too much on other parameters and data
   temp.R_L = C2K(14); units.temp.R_L = 'K'; label.temp.R_L = 'temperature';
   comment.R_L = 'based on Tinduff Hatchery observations, for the individuals of 5-6 cm length, nb of eggs for one spawning event';

   
 data.E_0 = 3e-3; units.E_0  = 'J';  label.E_0  = 'egg energy content'; bibkey.E_0  = 'RegnierBrisson';
  comment.E_0 = 'calculated with a 30% GSI .E0 = 1e-3J for Crassostrea gigas';  
  
data.GSI = 0.3; units.GSI  = '-';  label.GSI  = 'gonado-somatic index'; bibkey.GSI  = 'RegnierBrisson';
temp.GSI = C2K(14); units.temp.GSI = 'K'; label.temp.GSI = 'temperature';
 comment.GSI = 'calculated from 3-year individuals at spawning. We assume that they can spawn (at least) twice a year ';  
  


%data.Li  = 5.2;   units.Li  = 'cm';  label.Li  = 'ultimate shell height';   bibkey.Li  = 'RegnierBrisson';
%  comment.Li = 'mean obs collector at 3 years';
%   temp.Li = C2K(); units.temp.Li = 'K'; label.temp.Li = 'temperature';

% data.Ri  = 2*4000000/365; units.Ri  = '#/d'; label.Ri  = 'maximum reprod rate';     bibkey.Ri  = 'TinduffHatchery';   %test *2 for the Ri, 27/01/2023 -- impact of the kappa, but not too much on other parameters and data
%   temp.Ri = C2K(14); units.temp.Ri = 'K'; label.temp.Ri = 'temperature';
%   comment.Ri = 'based on Tinduff Hatchery observations, for the individuals of 5-6 cm length, nb of eggs for one spawning event';

 %data.Wdi = 1.6;    units.Wdi = 'g'; label.Wdi = 'ultimate dry weight';     bibkey.Wdi = 'RegnierBrisson';
 %  comment.Wdi = 'mean obs collector at 3 years';

% uni-variate data

%time-length larvae
tLTinduff = [ ... % time since spawning event (d). shell length (cm)
% 2   0.0105
5	0.01294
5	0.01326
5	0.01281
5	0.0129
5	0.01339
5	0.01251
5	0.01324
5	0.01319
5	0.01293
5	0.01316
5	0.01288
5	0.01295
5	0.01324
5	0.01298
5	0.01279
5	0.01305
5	0.01298
5	0.01263
5	0.01296
5	0.01289
5	0.0133
5	0.01331
5	0.01328
5	0.01286
5	0.01303
5	0.01316
5	0.01324
5	0.01304
5	0.01331
5	0.01301
7	0.01389
7	0.01453
7	0.01419
7	0.01419
7	0.01403
7	0.01414
7	0.01409
7	0.01381
7	0.01428
7	0.01399
7	0.01393
7	0.01456
7	0.0141
7	0.01419
7	0.01411
7	0.01423
7	0.0143
7	0.01436
7	0.01423
7	0.01443
7	0.01408
7	0.01385
7	0.01454
7	0.0139
7	0.01411
7	0.01449
7	0.01411
7	0.01423
7	0.01418
7	0.01449
9	0.01593
9	0.01565
9	0.01558
9	0.01568
9	0.01543
9	0.01556
9	0.01569
9	0.01529
9	0.01626
9	0.0159
9	0.01569
9	0.01551
9	0.01515
9	0.01558
9	0.01591
9	0.01528
9	0.01581
9	0.01544
9	0.01586
9	0.01508
9	0.01553
9	0.01561
9	0.01563
9	0.0159
9	0.01553
9	0.01578
9	0.01563
9	0.01546
9	0.01571
9	0.0158
11	0.01756
11	0.0174
11	0.01739
11	0.01714
11	0.01741
11	0.01736
11	0.01726
11	0.01771
11	0.01744
11	0.01664
11	0.01707
11	0.01747
11	0.01709
11	0.01763
11	0.0171
11	0.01741
11	0.01744
11	0.0166
11	0.01703
11	0.01699
11	0.01757
11	0.01731
11	0.01696
11	0.01737
11	0.01713
11	0.0172
11	0.01733
11	0.01673
11	0.01739
11	0.01734
13	0.0201
13	0.01933
13	0.01893
13	0.01836
13	0.01911
13	0.01823
13	0.01859
13	0.01865
13	0.0193
13	0.01788
13	0.01845
13	0.01844
13	0.01878
13	0.01803
13	0.01878
13	0.01883
13	0.01895
13	0.01868
13	0.01904
13	0.01875
13	0.01845
13	0.01875
13	0.0182
13	0.01935
13	0.01946
13	0.01925
13	0.0187
13	0.01853
13	0.01846
13	0.01851
15	0.01983
15	0.02014
15	0.02003
15	0.01994
15	0.0201
15	0.01999
15	0.0204
15	0.01974
15	0.02021
15	0.02023
15	0.02023
15	0.02006
15	0.01978
15	0.01973
15	0.01935
15	0.02019
15	0.0195
15	0.02013
15	0.01978
15	0.02015
15	0.02033
15	0.02009
15	0.02021
15	0.02004
15	0.01988
15	0.02035
15	0.01949
15	0.01939
15	0.02016
15	0.01998
17	0.02024
17	0.02052
17	0.01976
17	0.0204
17	0.02036
17	0.02008
17	0.02064
17	0.02022
17	0.0201
17	0.02033
17	0.0206
17	0.02038
17	0.02025
17	0.0197
17	0.02008
17	0.0202
17	0.02005
17	0.01973
17	0.0202
17	0.02063
17	0.0201
17	0.01957
17	0.0205
17	0.01977
17	0.0205
17	0.02025
17	0.0195
17	0.0199
17	0.0195
17	0.0207
% 18  0.0215  %aj     Lj 
30	0.044
30	0.0561
30	0.0506
30	0.051
30	0.0418
30	0.0522
30	0.0468
30	0.0435
30	0.0479
30	0.0607
30	0.047
30	0.0502
30	0.0543
30	0.041
30	0.0343
30	0.0578
30	0.0569
30	0.0604
30	0.0394
30	0.0659
30	0.0524
30	0.0617
30	0.0586
30	0.0396
30	0.0593
30	0.0561
30	0.0366
30	0.0372
30	0.0444
30	0.0574
35	0.0717
35	0.0797
35	0.0684
35	0.0709
35	0.078
35	0.0656
35	0.0768
35	0.0877
35	0.0637
35	0.0725
35	0.0656
35	0.0779
35	0.0778
35	0.0775
35	0.088
35	0.0799
35	0.0682
35	0.0649
35	0.0671
35	0.0618
35	0.0847
35	0.0793
35	0.0476
35	0.077
35	0.0574
35	0.0704
35	0.0677
35	0.0662
35	0.064
35	0.098
    ];

%%% differenciate larvae and juvenile from one data set
temp_l = tLTinduff(:,2);
% data.tLlarvae(:,1) = tLTinduff(tLTinduff(:,1) < data.aj); %only for values under the metamorphosis age
% data.tLlarvae(:,2) = temp_l(tLTinduff(:,1) < data.aj);
data.tLlarvae(:,1) = [data.ab ; tLTinduff(tLTinduff(:,1) < data.aj)]; %only for values under the metamorphosis age
data.tLlarvae(:,2) = [data.Lb ; temp_l(tLTinduff(:,1) < data.aj)];
units.tLlarvae   = {'d', 'cm'};  label.tLlarvae = {'time', 'shell length'};  %shell length, parallel to the hinge line, biggest length
temp.tLlarvae    = C2K(18);  units.temp.tLlarvae = 'K'; label.temp.tLlarvae = 'temperature'; %for the moment, 16°C to be the mean of temperature during nursery
bibkey.tLlarvae = 'TinduffHatchery';
comment.tLlarvae = 'Monitoring in the Tinduff Hatchery, controlled temperatures';

%%% Relation from Laure Regnier-Brisson data: Shell height = 1.0789 * shell length + 0.01905 (in cm), R² = 0.9841, N = 3208 individuals
% data.tLjuve(:,1) = tLTinduff(tLTinduff(:,1) >= data.aj); %only for values above the metamorphosis age
% data.tLjuve(:,2) = temp_l(tLTinduff(:,1) >= data.aj);
tLTinduff(:,2) = 1.1451 * tLTinduff(:,2); % 

data.tLjuve(:,1) = [data.aj ; tLTinduff(tLTinduff(:,1) >= data.aj)]; %only for values above the metamorphosis age
data.tLjuve(:,2) = [data.Lj ; temp_l(tLTinduff(:,1) >= data.aj)];
units.tLjuve   = {'d', 'cm'};  label.tLjuve = {'time', 'shell height'};  %shell height, perpendicular to the hinge line
temp.tLjuve    = C2K(18);  units.temp.tLjuve = 'K'; label.temp.tLjuve = 'temperature'; %for the moment, 16°C to be the mean of temperature during nursery
bibkey.tLjuve = 'TinduffHatchery';
comment.tLjuve = 'Monitoring in the Tinduff Hatchery, controlled temperatures';

%%% -----------------------------------------------------------------------
%%% Time - Shell height - Ash free total dry Weight - Ash free gonad dry weight
%%% Data from Laure Régnier-Brisson thesis, in Saint-Anne then Roscanvel
%%% The temperature is similar in the 2 places, so the real temperature of
%%% Saint-Anne will be considered.
%%% The temperature is defined with Fourier parameters determined in a
%%% previous step
%%% The translation shell length VS shell height (in mm) could be defined by:
            %%%     Shell length = 0.912 * Shell height - 1.309 for all size (N = 3208) and R² = 0.984
            %%%     Shell length = 0.843 * Shell height + 0.019 for height < 15 mm (N = 99) and R² = 0.972
%modified on 31/01/2023 - add the missing values with W as NaN and change
%the date to have the 0 at 24/05/2019 date of spawning
tLW19SA = [ %time since fertilisation (days), shell height (cm), ash free total dry weight (g), ash free gonad dry weight (g),
    % at 152 days, the organisms are still juveniles, but the shape is
    % considered the same for the whole data set
152	1.743	0.012	0
152	1.414	0.007	0
152	1.118	0.006	0
152	1.348	0.005	0
152	1.389	0.004	0
152	1.351	0.008	0
152	1.224	0.006	0
152	1.436	0.011	0
152	1.36	0.007	0
152	1.764	0.009	0
152	1.49	0.01	0
152	1.218	0.006	0
152	1.315	0.009	0
152	1.139	0.002	0
152	1.439	0.006	0
152	1.502	0.008	0
152	1.596	0.005	0
152	1.354	0.005	0
152	1.4	0.004	0
152	1.376	0.005	0
152	1.116	0.001	0
152	1.273	0.004	0
152	1.391	0.005	0
152	1.307	0.003	0
152	1.3	0.004	0
152	1.211	0.003	0
152	1.216	0.001	0
152	1.433	0.003	0
152	1.216	0.003	0
152	1.446	0.005	0
152	1.544	0.008	0
152	1.087	0.004	0
152	1.239	0.003	0
152	1.253	0.003	0
152	1.283	0.003	0
182	1.205	0.014	0
182	1.42	0.014	0
182	1.257	0.01	0
182	1.436	0.025	0
182	1.339	0.019	0
182	1.383	0.012	0
182	1.922	0.047	0
182	1.106	0.004	0
182	1.362	0.01	0
182	1.666	0.026	0
182	1.936	0.042	0
182	1.319	0.011	0
182	1.396	0.019	0
182	1.425	0.017	0
182	1.608	0.028	0
182	1.224	0.013	0
182	1.15	0.008	0
182	1.343	0.01	0
182	1.346	0.013	0
182	1.453	0.019	0
182	1.49	0.021	0
182	1.71	0.032	0
182	1.484	0.019	0
182	1.146	0.009	0
182	1.18	0.007	0
397	3.023	0.272	0
397	2.949	0.211	0
397	2.661	NaN	0
397	2.185	0.083	0
397	2.617	0.134	0
397	3.065	0.268	0
397	2.469	0.146	0
397	3.412	0.424	0
397	2.909	0.233	0
397	2.328	0.104	0
397	2.031	0.08	0
397	3.089	NaN	0
397	3.631	0.369	0
397	1.933	0.063	0
397	2.394	0.159	0
397	2.607	0.147	0
397	1.572	0.028	0
397	2.642	0.158	0
397	2.553	0.184	0
397	2.89	0.236	0
397	3.404	0.371	0
397	2.538	0.14	0
397	2.272	0.102	0
397	2.547	0.127	0
397	2.76	0.169	0
397	3.317	NaN	0
397	2.44	NaN	0
397	2.991	NaN	0
397	2.319	NaN	0
397	2.357	NaN	0
423	3.209	0.386	0.057
423	2.567	0.166	0.017
423	2.795	0.28	0.034
423	3.031	0.436	0.077
423	2.168	0.121	0.011
423	3.269	0.461	0.08
423	2.683	0.266	0.029
423	3.53	0.516	0.056
423	3.192	0.412	0.046
423	2.533	0.197	0.035
423	2.637	0.228	0.028
423	2.868	0.288	0.036
423	3.124	0.364	0.05
423	2.881	0.219	0.031
423	2.994	0.311	0.035
423	2.716	0.26	0.037
423	3.484	0.522	0.076
423	3.103	0.419	0.077
423	2.954	0.299	0.033
423	3.016	0.34	0.066
423	3.174	0.405	0.04
423	2.9	0.355	0.04
423	2.254	NaN	NaN
423	2.747	0.254	0.036
423	3.122	0.393	0.053
423	3.188	0.334	0.044
423	2.988	0.33	0.033
423	3.112	0.486	0.001
423	3.016	0.4	0.056
423	3.412	0.508	0.1
451	3.057	NaN	0.064
451	2.26	NaN	0.013
451	3.417	NaN	0.087
451	3.128	NaN	0.045
451	3.009	NaN	0.05
451	3.186	0.418	0.053
451	3.34	0.415	0.056
451	3.392	0.485	0.07
451	2.881	0.27	0.022
451	3.268	0.347	0.049
451	3.126	0.39	0.05
451	3.647	0.475	0.063
451	2.86	0.259	0.025
451	3.026	0.34	0.042
451	3.875	0.726	0.089
451	3.261	0.414	0.07
451	3.098	0.366	0.058
451	3.224	0.442	0.095
451	2.905	0.333	0.048
451	2.544	0.083	0.005
451	2.92	0.243	0.019
451	3.358	0.473	0.066
451	3.372	0.359	0.047
451	3.164	0.38	0.059
451	2.657	0.208	0.036
451	3.513	NaN	0.073
451	3.095	NaN	0.041
451	3.94	NaN	0.075
451	2.944	NaN	0.068
451	3.776	NaN	0.134
479	3.676	NaN	0.054
479	3.768	NaN	0.038
479	3.99	NaN	0.075
479	3.644	NaN	0.064
479	3.512	NaN	0.036
479	3.65	0.443	0.032
479	3.577	0.461	0.028
479	3.32	0.392	0.011
479	3.524	0.449	0.0090
479	4.198	0.765	0.039
479	3.453	0.418	0.037
479	4.3	0.822	0.017
479	3.574	0.62	0.048
479	3.38	0.421	0.022
479	3.416	0.554	0.046
479	2.826	0.333	0.013
479	4.629	0.905	0.072
479	4.093	0.861	0.034
479	4.353	0.869	0.033
479	3.795	0.528	0.05
479	4.007	0.706	0.03
479	3.381	0.422	0.024
479	4.271	0.903	0.044
479	3.885	0.357	0.02
479	3.437	0.513	0.013
479	3.752	NaN	0.024
479	3.786	NaN	0.023
479	3.775	NaN	0.084
479	3.447	NaN	0.047
479	4.24	NaN	0.066
508	3.592	NaN	0.012
508	3.407	NaN	0.0090
508	3.528	NaN	0.012
508	2.799	NaN	0.011
508	2.914	NaN	0.0060
508	3.846	0.724	0.033
508	3.62	0.506	0.028
508	3.093	0.322	0.012
508	3.95	0.714	0.014
508	4.158	0.729	0.012
508	4.332	0.749	0.024
508	4.273	0.835	0.021
508	3.759	0.51	0.013
508	4.239	0.928	0.021
508	3.811	0.574	0.015
508	3.138	0.432	0.013
508	3.866	0.653	0.017
508	3.937	0.685	0.017
508	3.75	0.597	0.014
508	4.123	0.681	0.016
508	3.551	0.542	0.011
508	3.545	0.315	0.015
508	3.134	0.331	0.0090
508	3.132	0.355	0.011
508	3.033	0.266	0.005
508	3.532	NaN	0.015
508	3.369	NaN	0.0090
508	4.205	NaN	0.045
508	3.944	NaN	0.022
508	3.386	NaN	0.013
536	3.946	NaN	0.012
536	3.441	NaN	0.015
536	4.655	NaN	0.014
536	3.887	NaN	0.014
536	3.857	NaN	0.01
536	4.109	0.528	0.018
536	3.12	0.22	0.005
536	4.088	0.357	0.0070
536	4.294	0.708	0.012
536	4.222	0.497	0.0080
536	4.036	0.411	0.015
536	3.463	0.574	0.019
536	3.86	0.852	0.02
536	3.635	0.374	0.013
536	3.217	0.798	0.018
536	3.811	0.816	0.018
536	3.246	0.841	0.01
536	4.436	0.794	0.016
536	4.582	0.688	0.025
536	3.709	0.301	0.0090
536	4.098	0.722	0.014
536	3.831	0.657	0.016
536	3.703	0.882	0.017
536	4.136	0.777	0.0080
536	3.805	0.818	0.017
536	4.481	NaN	0.02
536	3.472	NaN	0.0080
536	3.27	NaN	0.017
536	4.431	NaN	0.026
536	3.752	NaN	0.024
593	3.853	NaN	0.0090
593	4.517	NaN	0.013
593	3.879	NaN	0.01
593	3.999	NaN	0.011
593	4.471	NaN	0.011
593	4.599	0.697	0.02
593	4.596	0.821	0.021
593	4.147	0.761	0.019
593	4.6	0.735	0.014
593	4.362	0.653	0.005
593	4.371	0.854	0.012
593	4.61	0.804	0.024
593	3.841	0.476	0.011
593	4.345	0.701	0.015
593	4.365	0.841	0.015
593	3.652	0.409	0.0060
593	4.545	0.748	0.012
593	4.176	0.518	0.011
593	4.009	0.577	0.015
593	3.023	0.259	0.0080
593	4.153	0.493	0.013
593	4.929	0.979	0.02
593	4.561	0.761	0.019
593	4.498	0.579	0.0060
593	4.694	0.858	0.014
593	4.266	NaN	0.014
593	3.199	NaN	0.005
593	3.509	NaN	0.012
593	4.652	NaN	0.018
593	4.597	NaN	0.014
654	4.131	NaN	0.032
654	4.562	NaN	0.056
654	3.78	NaN	0.015
654	4.338	NaN	0.018
654	4.464	NaN	0.027
654	4.447	0.984	0.03
654	4.498	0.883	0.026
654	3.934	0.393	0.014
654	4.814	1.025	0.21
654	4.238	0.662	0.015
654	3.918	0.569	0.02
654	4.356	0.91	0.049
654	4.695	0.822	0.017
654	3.979	0.52	0.021
654	5.123	1.387	0.097
654	4.114	0.488	0.018
654	3.837	0.318	0.011
654	4.251	0.715	0.019
654	4.463	0.674	0.016
654	4.849	0.848	0.025
654	4.326	0.614	0.018
654	4.142	0.595	0.028
654	5.167	1.153	0.033
654	4.694	0.731	0.027
654	4.348	0.862	0.077
654	4.118	NaN	0.056
654	3.806	NaN	0.015
654	4.488	NaN	0.019
654	4.12	NaN	0.043
654	4.259	NaN	0.018
1189	4.885	0.748829	0.275445
1189	4.186	0.665021	0.208701
1189	4.634	0.994237	0.295856
1189	5.398	1.568337	0.468056
1189	5.246	1.604477	0.393071
1189	5.183	1.305273	0.386319
1189	4.522	1.159030	0.250403
1189	5.818	1.714369	0.506603
1189	5.512	1.894193	0.554741
1189	4.378	0.806979	0.23225
1189	5.313	1.273014	0.418932
1189	5.573	1.833609	0.504891
1189	5.475	1.936551	0.517255
1189	5.26	1.307569	0.38458
1189	5.246	1.621176	0.443907
1189	5.007	1.035512	0.313676
1189	5.113	1.201195	0.380353
1189	5.986	1.894551	0.524977
1189	4.794	0.868941	0.282705
1189	5.518	2.013026	0.47769
1189	5.838	2.312711	0.575377
1189	5.702	1.882206	0.550472
1189	5.059	1.599895	0.439025
1189	5.327	1.613184	0.447415
1189	5.721	1.823543	0.549468
1189	5.447	2.076962	0.478787
1189	5.354	1.609162	0.452883
1189	5.109	1.7056      0.355435
1189	5.956	2.027093	0.554408
1189	5.413	1.837877	0.495532
];
%%% to delete NaN values, in case the matlab version does not allow to deal
%%% with them
tLW19SA(any(isnan(tLW19SA), 2), :) = [];

% data.tL19SA(:,1) = tLW19SA(:,1); data.tL19SA(:,2) = tLW19SA(:,2); %creat time-length data
data.tL19SA(:,1) = [data.aj ; tLW19SA(:,1)]; data.tL19SA(:,2) = [data.Lj ; tLW19SA(:,2)]; %creat time-length data with th initial data as the metamorphosis data
units.tL19SA   = {'d', 'cm'};  label.tL19SA = {'time', 'shell height', 'N19 at Ste Anne'};  
temp.tL19SA    = [  365         13.712             % Fourier parameters to be used in the derivatives function, from the excel file in Enviro/
                    -0.76009	-3.7422             % Fourier parameters estimated for Saint-Anne, because the Roscanvel temperatures are similar starts on the 1st June 2019
                    -0.16068	0.031801
                    -0.056534	-0.073573];  units.temp.tL19SA = 'C'; label.temp.tL19SA = 'Fourier parameters for temperature cycle for 19SA data';
bibkey.tL19SA = 'RegnierBrisson';

% data.tW19SA(:,1) = tLW19SA(:,1);  data.tW19SA(:,2) = tLW19SA(:,3);
data.tW19SA(:,1) = [tLW19SA(1,1) ; tLW19SA(:,1)];  data.tW19SA(:,2) = [tLW19SA(1,3) ; tLW19SA(:,3)]; %repeat the first row to have the same dimension than tL
units.tW19SA   = {'d', 'g'};  label.tW19SA = {'time', 'ash free dry weight', 'N19 at Ste Anne'};  
temp.tW19SA    = [  365         13.712             % Fourier parameters to be used in the derivatives function, from the excel file in Enviro/
                    -0.76009	-3.7422             % Fourier parameters estimated for Saint-Anne, because the Roscanvel temperatures are similar
                    -0.16068	0.031801
                    -0.056534	-0.073573];  units.temp.tW19SA = 'C'; label.temp.tW19SA = 'Fourier parameters for temperature cycle for 19SA data';
bibkey.tW19SA = 'RegnierBrisson';


data.LW19SA(:,1) = tLW19SA(:,2); data.LW19SA(:,2) = tLW19SA(:,3) - tLW19SA(:,4); %creat time-length data
units.LW19SA   = {'cm', 'g'};  label.LW19SA = {'shell height','ash free dry weight without gonad dry weight', 'N19 at Ste Anne'};  
temp.LW19SA    = C2K(14) ;  units.temp.LW19SA = 'C'; label.temp.LW19SA = 'temperature';
bibkey.LW19SA = 'RegnierBrisson';
comment.LW19SA = 'temperature between 17 and 9 C - Mean temperature (monitoring)';




%% set weights for all real data
weights = setweights(data, []);

weights.tLlarvae        = 5 * weights.tLlarvae;
weights.tLjuve          = 0 * weights.tLjuve;
weights.tL19SA        = 5 * weights.tL19SA; 
% weights.tL19SA(1:60)    = 1 * weights.tL19SA(1:60); 
% weights.tL19SA(61:end)  = 5 * weights.tL19SA(61:end); 
weights.tW19SA(1)       = 0 * weights.tW19SA(1); % for the first value that is repeated just to have the right dimension
weights.tW19SA(2:end)   = 4 * weights.tW19SA(2:end); 
weights.LW19SA(1:60)    = 2 * weights.LW19SA(1:60); %can add weight for the first values that are never fitted corectly
weights.LW19SA(61:end)  = 2 * weights.LW19SA(61:end); 

weights.Lb = 2 * weights.Lb;
weights.Lj = 1 * weights.Lj;
weights.Lp = 2 * weights.Lp; 
%weights.Li = 2 * weights.Li;

weights.ab = 2 * weights.ab;
weights.aj = 1 * weights.aj;
weights.tp = 0 * weights.tp;
weights.am = 1 * weights.am;

%weights.Ri = 1 * weights.Ri;
weights.R_L = 1 * weights.R_L;
weights.E_0 = 10 * weights.E_0;


%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);
data.psd.kap = 0.88;    units.psd.kap = '-';        label.psd.kap = 'allocation fraction to soma for Pectinids';
weights.psd.kap = 80 * weights.psd.kap;
% weights.psd.v = 5 * weights.psd.v;
%weights.psd.p_M = 10 * weights.psd.p_M;
% data.psd.del_M = 0.3; units.psd.k  = '-'; label.psd.k  = 'shape coefficient'; 
% weights.psd.del_M = 5 * weights.psd.del_M;

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;

% %% Group plots
% set1 = {'tL19SA'}; subtitle1 = {'Data for SA N19'};
% metaData.grp.sets = {set1};
% metaData.grp.subtitle = {subtitle1};


%% Discussion points
D1 = 'Temperature (C) of tL data changes as: T(t)=mean+amplitude*sin(2*pi*(t+250)/365)';
metaData.discussion = struct('D1',D1);

%% Facts
F1 = 'Incidental hermaphroditic';
metaData.bibkey.F1 = 'Wiki';
F2 = 'Sports jet propulsion swimming';
metaData.bibkey.F2 = 'Wiki';
metaData.facts = struct('F1',F1, 'F2',F2);

%% Links
metaData.links.id_CoL = '6RKGL'; % Cat of Life
metaData.links.id_ITIS = '79628'; % ITIS
metaData.links.id_EoL = '46468063'; % Ency of Life
metaData.links.id_Wiki = 'Mimachlamys'; % Wikipedia
metaData.links.id_ADW = ''; % ADW
metaData.links.id_Taxo = '3967485'; % Taxonomicon
metaData.links.id_WoRMS = '236719'; % WoRMS
metaData.links.id_molluscabase = '236719'; % molluscabase


%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
'howpublished = {\url{http://en.wikipedia.org/wiki/Chlamys_varia}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
'author = {Kooijman. S.A.L.M.}. ' ...
'year = {2010}. ' ...
'title  = {Dynamic Energy Budget theory for metabolic organisation}. ' ...
'publisher = {Cambridge Univ. Press. Cambridge}. ' ...
'pages = {Table 4.2 (page 150). 8.1 (page 300)}. ' ...
'howpublished = {\url{../../../bib/Kooy2010.html}}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'marlin'; type = 'Misc'; bib = ...
'howpublished = {\urlhttp://www.marlin.ac.uk/biotic/browse.php?sp=6255}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'sealifebase'; type = 'Misc'; bib = ...
'howpublished = {\url{http://www.sealifebase.org/summary/Mimachlamys-varia.html}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'TinduffHatchery'; type = 'Misc'; bib = ...
'TinduffHatchery';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'RegnierBrisson'; type = 'Misc'; bib = ...
'RegnierBrisson';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
