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

function [icaS,icaA,correlations_dict] = ICAdim2_noPCA(X, C);

CorrThresh=0.9;     % reduce this to make split-half reproducibility less conservative default=0.9
RepeatsPerJ=10;     % number of random re-run repeats for each PCA dimensionality
g='tanh';           % fastica 'g' option - typically 'tanh' or 'pow3'.  default=tanh
epsilon=1e-13;      % fastica convergence threshold. default=1e-13
DoPass2=0;          % do we run the final optimisation at the chosen dimensionalities
ShowRandomTries=0;  % how many of the random reruns to show debugging plot for. Set to 0 for none.
vrb = 'off';        % indicate if the fastica function should be verbose or not
itr = 10000;         % number of iterations. default = 3000
aprc='symm';        % The decorrelation approach used: symmetric ('symm'), or deflation ('defl'). default='symm'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if length(Jvalues==1) & Jvalues<0
%   J=-Jvalues;
%   [pcaU,pcaS,pcaV]=nets_svds(X,J);
%   [icaS,icaA,~]=fastica(pcaV','approach',aprc,'g',g,'epsilon',epsilon,'maxNumIterations',itr, 'verbose', vrb);
%   return  % I'm so so so sorry
% end

%%% first-pass: compare a range of dimensionalities
% [pcaU,pcaS,pcaV]=nets_svds(X,max(Jvalues));
NStrongCorrs=[];
NStrongCorrs_sum = [];
correlations_dict = containers.Map;
% for J=Jvalues
J = min(size(X));
grot2=[];
for JJ=1:RepeatsPerJ
if C==1
  grot=randn(size(X,1),1)>0;  X1 = X(grot==0,:);  X2 = X(grot==1,:);
else
  grot=randn(size(X,2),1)>0;  X1 = X(:,grot==0);  X2 = X(:,grot==1);
end
[icaS,icaA,~]   = fastica(X,            'approach',aprc,'g',g,'epsilon',epsilon,'maxNumIterations',itr, 'verbose', vrb,'lastEig',J);
[icaS1,icaA1,~] = fastica(X1,       'approach',aprc,'g',g,'epsilon',epsilon,'maxNumIterations',itr, 'verbose', vrb,'lastEig',J);
[icaS2,icaA2,~] = fastica(X2,       'approach',aprc,'g',g,'epsilon',epsilon,'maxNumIterations',itr, 'verbose', vrb,'lastEig',J);
if C==1
%       if (size(icaS,1) ~= J) || (size(icaS1,1) ~= J) || (size(icaS2,1) ~= J)
%           break
%       end
  [grotA,grotB,grotC,grotD,grotE,grotF] = ssorder2(icaS',icaA',icaS1',icaA1',icaS2',icaA2');
  icaS=grotA';icaA=grotB';icaS1=grotC';icaA1=grotD';icaS2=grotE';icaA2=grotF';
  grot=corr([icaS' icaS1' icaS2']); % match on ICA sources
else
%       icaA=pcaU(:,1:J)*pcaS(1:J,1:J)*icaA; icaA1=pcaU1*pcaS1*icaA1; icaA2=pcaU2*pcaS2*icaA2;
  [icaA,icaS,icaA1,icaS1,icaA2,icaS2] = ssorder2(icaA,icaS,icaA1,icaS1,icaA2,icaS2);
  grot=corr([icaA icaA1 icaA2]); % match on ICA mixing matrix, projected back to input space for compatibility across folds

end
grot2=[grot2;diag(grot(J+1:2*J,2*J+1:end))'];  % focus on similarity between the two split-halves
%     if JJ <= ShowRandomTries   % show individual sub-runs
%       JJJ=min(J,34);
%       figure; n=ceil(sqrt(JJJ+2));
%         for i=1:JJJ
%           subplot(n,n,i);
%           if C==1
%             plot([icaS(i,:)' icaS1(i,:)' icaS2(i,:)']);
%           else
%             plot([icaA(:,i) icaA1(:,i) icaA2(:,i)]);
%           end
%         end
%       subplot(n,n,JJJ+1); imagesc(grot,[-1 1]); subplot(n,n,JJJ+2); plot(diag(grot(J+1:2*J,2*J+1:end))');
%     end
end;
figure; plot(grot2'+randn(size(grot2'))*0.001); hold on; plot(mean(grot2)','k','LineWidth',2); title(sprintf('J=%d',J)); drawnow;
correlations_dict(int2str(J)) = grot2;
NStrongCorrs=[NStrongCorrs sum(mean(grot2)>CorrThresh)];
grot2=mean(grot2);  NStrongCorrs_sum=[NStrongCorrs_sum sum(grot2(grot2>CorrThresh))];
% end;
% figure; plot(J,NStrongCorrs); title('Number of Strong Correlations'); drawnow;
% figure; plot(J,NStrongCorrs_sum); title('Sum of Strong Correlations');drawnow;
NStrongCorrs

% if DoPass2
%   %%% second-pass: find best ICA run for chosen dimensionality
%   grot=find(NStrongCorrs==max(NStrongCorrs)); grot=grot(1)
%   J=Jvalues(grot);               % PCA dimensionality input to ICA
%   JJJ=ceil(NStrongCorrs(grot));  % how many ICA components we will keep
%   grot2=[]; grot3=-1;
%   for JJ=1:RepeatsPerJ*3
%       if C==1
%         grot=randn(size(X,1),1)>0;  [pcaU1,pcaS1,pcaV1]=nets_svds(X(grot==0,:),J);  [pcaU2,pcaS2,pcaV2]=nets_svds(X(grot==1,:),J);
%       else
%         grot=randn(size(X,2),1)>0;  [pcaU1,pcaS1,pcaV1]=nets_svds(X(:,grot==0),J);  [pcaU2,pcaS2,pcaV2]=nets_svds(X(:,grot==1),J);
%       end
%       [icaS,icaA,~]   = fastica(pcaV(:,1:J)', 'approach',aprc,'g',g,'epsilon',epsilon,'maxNumIterations',itr, 'verbose', vrb);
%       [icaS1,icaA1,~] = fastica(pcaV1',       'approach',aprc,'g',g,'epsilon',epsilon,'maxNumIterations',itr, 'verbose', vrb);
%       [icaS2,icaA2,~] = fastica(pcaV2',       'approach',aprc,'g',g,'epsilon',epsilon,'maxNumIterations',itr, 'verbose', vrb);
%       if C==1
%         [grotA,grotB,grotC,grotD,grotE,grotF] = ssorder2(icaS',icaA',icaS1',icaA1',icaS2',icaA2');
%         icaS=grotA';icaA=grotB';icaS1=grotC';icaA1=grotD';icaS2=grotE';icaA2=grotF';
%         grot=corr([icaS' icaS1' icaS2']);
%       else
%         icaA=pcaU(:,1:J)*pcaS(1:J,1:J)*icaA; icaA1=pcaU1*pcaS1*icaA1; icaA2=pcaU2*pcaS2*icaA2;
%         [icaA,icaS,icaA1,icaS1,icaA2,icaS2] = ssorder2(icaA,icaS,icaA1,icaS1,icaA2,icaS2);
%         grot=corr([icaA icaA1 icaA2]); 
%       end
%       grot2=[grot2 mean(diag(grot(J+1:J+JJJ,2*J+1:2*J+JJJ)))];
%       if grot2(end)>grot3
%         grot3=grot2(end);  icaAout=icaA(:,1:JJJ);  icaSout=icaS(1:JJJ,:);      
%         if C==2
%           icaAout=inv(pcaS(1:J,1:J))*pcaU(:,1:J)'*icaAout; % revert back to ICA-space mixing matrix
%         end
%       end
%   end;
%   icaA=icaAout; icaS=icaSout;
% end;

J

