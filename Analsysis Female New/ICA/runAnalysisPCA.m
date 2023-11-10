addpath FSLNets FastICA_25 FACS labelpoints

load('M_deltas_deconf.mat')
M = X_deconf;

X = M;
clear label

X = nets_normalise(X);

X = X';
disp(size(X));
Jvalues = [2:57];
[pcaU,pcaS,pcaV]=nets_svds(X, max(Jvalues));

figure;
plot(diag(pcaS));
xlabel('Singurlar Value Number/Rank');
ylabel('Singural Value Intensity');
grid('on');
grid minor;

percVariance = cumsum(diag(pcaS).^2)*100 / sum(diag(pcaS).^2);

figure;
plot(percVariance);
xlabel('Singurlar Value Number/Rank');
ylabel('Percentage (%) Variance');
grid('on');
grid minor;

figure;
scatter(pcaU(:,1), pcaU(:,2))
labels = {'T1_nonlinear', 'T1_linear', 'jacobian', 'vbm', 'T2_nonlinear', 'T2_lesions', 'swi','rsfmri_0', 'rsfmri_1', 'rsfmri_2', 'rsfmri_3', 'rsfmri_4', 'rsfmri_5', 'rsfmri_6', 'rsfmri_7', 'rsfmri_8','rsfmri_9', 'rsfmri_10', 'rsfmri_11', 'rsfmri_12', 'rsfmri_13', 'rsfmri_14', 'rsfmri_15','rsfmri_16', 'rsfmri_17', 'rsfmri_18', 'rsfmri_19', 'rsfmri_20', 'rsfmri_21', 'rsfmri_22','rsfmri_23', 'rsfmri_24', 'tfmri_1', 'tfmri_2', 'tfmri_5', 'tfmri_c_1', 'tfmri_c_2', 'tfmri_c_5','tracts', 'tbss_FA_s', 'tbss_ICVF_s', 'tbss_ISOVF_s', 'tbss_L1_s', 'tbss_L2_s', 'tbss_L3_s', 'tbss_MD_s', 'tbss_MO_s', 'tbss_OD_s', 'tbss_FA', 'tbss_ICVF', 'tbss_ISOVF', 'tbss_L1', 'tbss_L2','tbss_L3', 'tbss_MD', 'tbss_MO', 'tbss_OD'};
h = labelpoints(pcaU(:,1), pcaU(:,2), labels);
grid on;
grid minor;
xlabel(['PCA Component 1 (' , num2str(percVariance(1)) , '%)'])
ylabel(['PCA Component 2 (' , num2str(percVariance(2) - percVariance(1)) , '%)'])

figure;
scatter(pcaV(:,1), pcaV(:,2))
grid on;
grid minor;
xlabel(['PCA Component 1 (' , num2str(percVariance(1)) , '%)'])
ylabel(['PCA Component 2 (' , num2str(percVariance(2) - percVariance(1)) , '%)'])

% % % % A = U S Vh'
% % % % A = m-by-n
% % % % U = m-by-p
% % % % S = p-by-p
% % % % Vh = p-by-n


X = X';
disp(size(X));
Jvalues = [2:50];
[pcaU,pcaS,pcaV]=nets_svds(X, max(Jvalues));

figure;
plot(diag(pcaS));
xlabel('Singurlar Value Number/Rank');
ylabel('Singural Value Intensity');
grid('on');
grid minor;
% 
percVariance = cumsum(diag(pcaS).^2)*100 / sum(diag(pcaS).^2);

figure;
plot(percVariance);
xlabel('Singurlar Value Number/Rank');
ylabel('Percentage (%) Variance');
grid('on');
grid minor;

figure;
scatter(pcaU(:,1), pcaU(:,2))
grid on;
grid minor;
xlabel(['PCA Component 1 (' , num2str(percVariance(1)) , '%)'])
ylabel(['PCA Component 2 (' , num2str(percVariance(2) - percVariance(1)) , '%)'])

figure;
scatter(pcaV(:,1), pcaV(:,2))
labels = {'T1_nonlinear', 'T1_linear', 'jacobian', 'vbm', 'T2_nonlinear', 'T2_lesions', 'swi','rsfmri_0', 'rsfmri_1', 'rsfmri_2', 'rsfmri_3', 'rsfmri_4', 'rsfmri_5', 'rsfmri_6', 'rsfmri_7', 'rsfmri_8','rsfmri_9', 'rsfmri_10', 'rsfmri_11', 'rsfmri_12', 'rsfmri_13', 'rsfmri_14', 'rsfmri_15','rsfmri_16', 'rsfmri_17', 'rsfmri_18', 'rsfmri_19', 'rsfmri_20', 'rsfmri_21', 'rsfmri_22','rsfmri_23', 'rsfmri_24', 'tfmri_1', 'tfmri_2', 'tfmri_5', 'tfmri_c_1', 'tfmri_c_2', 'tfmri_c_5','tracts', 'tbss_FA_s', 'tbss_ICVF_s', 'tbss_ISOVF_s', 'tbss_L1_s', 'tbss_L2_s', 'tbss_L3_s', 'tbss_MD_s', 'tbss_MO_s', 'tbss_OD_s', 'tbss_FA', 'tbss_ICVF', 'tbss_ISOVF', 'tbss_L1', 'tbss_L2','tbss_L3', 'tbss_MD', 'tbss_MO', 'tbss_OD'};
h = labelpoints(pcaV(:,1), pcaV(:,2), labels);
grid on;
grid minor;
xlabel(['PCA Component 1 (' , num2str(percVariance(1)) , '%)'])
ylabel(['PCA Component 2 (' , num2str(percVariance(2) - percVariance(1)) , '%)'])

tempdir = 'pca_variance_weights';
mkdir pca_variance_weights;
FolderName = tempdir;   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  savefig(fullfile(FolderName, [FigName '.fig']));
end
