%
% dimica - fastica wrapper, with cross-validated optimisation of both PCA and ICA dimensionality.
% Stephen Smith, FMRIB, Oxford University, September 2020.
% Maximises number of robustly-found ICA components.
% Note there is no fully left-out data (things are optimised for this data using this data).
%
% [icaS,icaA,pcaU,pcaS,pcaV] = ICAdim2(X,Jvalues,C);
%
% X is data matrix; ICA will find *horizontal* sources
%
% Jvalues is a horizontal vector of dimension values to try, e.g.:  10 or [10:5:100]
% if Jvalues is a single negative value, just run ICA once for this value.
%
% C is cross-validation direction 1=vertical 2=horizontal
% E.g., if you have X=subjectsXfeatures matrix:
%  - for feature-ICA:  ICAdim2(X,....,1);
%  - for subject-ICA:  ICAdim2(X',....,2);
%

function [icaS,icaA,pcaU,pcaS,pcaV] = ICAdim2_Pass2(X, Jvalues, correlations_dict,C);

CorrThresh=0.7;     % reduce this to make split-half reproducibility less conservative default=0.9
RepeatsPerJ=100;     % number of random re-run repeats for each PCA dimensionality
g='tanh';           % fastica 'g' option - typically 'tanh' or 'pow3'.  default=tanh
epsilon=1e-13;      % fastica convergence threshold. default=1e-13
DoPass2=1;          % do we run the final optimisation at the chosen dimensionalities
ShowRandomTries=0;  % how many of the random reruns to show debugging plot for. Set to 0 for none.
vrb = 'off';        % indicate if the fastica function should be verbose or not
itr = 10000;         % number of iterations. default = 3000
aprc='symm';        % The decorrelation approach used: symmetric ('symm'), or deflation ('defl'). default='symm'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pcaU,pcaS,pcaV]=nets_svds(X,max(Jvalues));

% load('PcaIca_subject_deltas_transpose_steve_norm.mat');
vs = values(correlations_dict); ks=keys(correlations_dict);
vs1 = []; ks1 = [];
idx=1;
for v=vs
ks1 = [ks1  str2double(cell2mat(ks(idx)))];
idx=idx+1;
vs1 = [vs1 sum(mean(cell2mat(v))>CorrThresh)];
end
[ks1_sorted, a_order] = sort(ks1);
vs1_sorted = vs1(a_order);
figure; plot(ks1_sorted,vs1_sorted); title(sprintf('Number of Strong correlations at thr=%d',CorrThresh)); drawnow;

grot=find(vs1_sorted==max(vs1_sorted)); grot=grot(1)
J=ks1_sorted(grot);               % PCA dimensionality input to ICA
JJJ=ceil(vs1_sorted(grot));  % how many ICA components we will keep

grot2=[]; grot3=-1;
for JJ=1:RepeatsPerJ*3
  if C==1
    grot=randn(size(X,1),1)>0;  [pcaU1,pcaS1,pcaV1]=nets_svds(X(grot==0,:),J);  [pcaU2,pcaS2,pcaV2]=nets_svds(X(grot==1,:),J);
  else
    grot=randn(size(X,2),1)>0;  [pcaU1,pcaS1,pcaV1]=nets_svds(X(:,grot==0),J);  [pcaU2,pcaS2,pcaV2]=nets_svds(X(:,grot==1),J);
  end
  [icaS,icaA,~]   = fastica(pcaV(:,1:J)', 'approach',aprc,'g',g,'epsilon',epsilon,'maxNumIterations',itr, 'verbose', vrb);
  [icaS1,icaA1,~] = fastica(pcaV1',       'approach',aprc,'g',g,'epsilon',epsilon,'maxNumIterations',itr, 'verbose', vrb);
  [icaS2,icaA2,~] = fastica(pcaV2',       'approach',aprc,'g',g,'epsilon',epsilon,'maxNumIterations',itr, 'verbose', vrb);
  if C==1
    [grotA,grotB,grotC,grotD,grotE,grotF] = ssorder2(icaS',icaA',icaS1',icaA1',icaS2',icaA2');
    icaS=grotA';icaA=grotB';icaS1=grotC';icaA1=grotD';icaS2=grotE';icaA2=grotF';
    grot=corr([icaS' icaS1' icaS2']);
  else
    icaA=pcaU(:,1:J)*pcaS(1:J,1:J)*icaA; icaA1=pcaU1*pcaS1*icaA1; icaA2=pcaU2*pcaS2*icaA2;
    [icaA,icaS,icaA1,icaS1,icaA2,icaS2] = ssorder2(icaA,icaS,icaA1,icaS1,icaA2,icaS2);
    grot=corr([icaA icaA1 icaA2]); 
  end
  grot2=[grot2 mean(diag(grot(J+1:J+JJJ,2*J+1:2*J+JJJ)))];
  if grot2(end)>grot3
    grot3=grot2(end);  icaAout=icaA(:,1:JJJ);  icaSout=icaS(1:JJJ,:);      
    if C==2
      icaAout=inv(pcaS(1:J,1:J))*pcaU(:,1:J)'*icaAout; % revert back to ICA-space mixing matrix
    end
  end
end;
icaA=icaAout; icaS=icaSout;

J, JJJ
