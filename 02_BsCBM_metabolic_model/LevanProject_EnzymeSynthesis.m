% This script was made to run all metabolic model for CASB-Levan Project.
% Created by Muhammad Naufal Hakim, Bioinformatics and Systems Biology KMUTT

% Don't forget to initiate COBRA Toolbox prior to run the script.

% To use this script, please look at modified excel file "iYO844_levansucrase_alpha_transition.xlsx".
% The sucrose input and alpha value is really depending on dynamic model part.

% Define the model from excel file
iYO844_enzyme_breakloop = xls2model...
    ('iYO844_levansucrase_alpha_transition.xlsx');

% FLUX BALANCE ANALYSIS USING MINIMAL TOTAL FLUX

% Flux simulation by using minimization total flux routine (Holzhutter et al, 2004; Lewis et al, 2010)
[MinimizedFlux, modelIrrev]= ...
    minimizeModelFlux(iYO844_enzyme_breakloop, 'min');

% Print flux simulation into designated matrix
reactions_check = modelIrrev.rxns;
flux_check = round(MinimizedFlux.x,10);
T_check = table(reactions_check,flux_check,...
    'RowNames',modelIrrev.rxns);
	
% Print the matrix into CSV file
writetable(T_check,...
    'Bsub_iYO844_breakloop_sucrose_uptake15.316_flux_2E-5_levansucrase_nh4_LB4.7.csv');

