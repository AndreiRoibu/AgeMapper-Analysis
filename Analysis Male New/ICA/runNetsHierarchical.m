% load('correlation_matrix_debiased_deconfounded.mat')
load('correlation_matrix_debiased_deconfounded_testA.mat')
addpath FSLNets

[dpRSN,yyRSN] = nets_hierarchy_andrei_mod(correlation_matrix_debiased_deconfounded);