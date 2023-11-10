
% nets_unconfound3(y,conf,keep) 
% nets_unconfound3(y,conf,keep,-1) 
%
% regresses conf out of y, handling missing data, and leaving the unique space of "keep" in
%
% data, confounds and output are all demeaned unless the "-1" option is included

function yd=nets_unconfound3(y,conf,keep,varargin)

DEMEAN=1;
if nargin==4
  DEMEAN=0;
end

if DEMEAN
  y = nets_demean(y);
  conf = nets_demean(conf);
  keep = nets_demean(keep);
end

if sum(isnan(conf(:))) == 0   % if there's no missing data in conf, can use faster code

  yd=y;
  yd(isnan(y))=0; % this is allowed because setting NaNs to zero (in y) doesn't change the maths...much....
  %conf=nets_svds(conf,rank(conf)); beta=pinv(conf)*yd; beta(abs(beta)<1e-10)=0; yd=y-conf*beta; % re-inherit nans
  conf=nets_svds(conf,rank(conf)); keep=nets_svds(keep,rank(keep));
  beta=pinv([conf keep])*yd; beta(abs(beta)<1e-10)=0; yd=y-conf*beta(1:size(conf,2),:); % re-inherit nans

  if DEMEAN
    yd = nets_demean( yd );
  end

else

  yd = (y*0)/0; % set to all NaN because we are not going to necessarily write into all elements below

  for i=1:size(y,2)
    grot=~isnan(sum([y(:,i) conf keep],2));
    grotconf = conf(grot,:); grotkeep = keep(grot,:);
    if DEMEAN
      grotconf=nets_demean(grotconf);
      grotkeep=nets_demean(grotkeep);
    end

    %%%grotconf=nets_svds(grotconf,rank(grotconf)); beta=pinv(grotconf)*y(grot,i); beta(abs(beta)<1e-10)=0; yd(grot,i)=y(grot,i)-grotconf*beta;
    grotconf=nets_svds(grotconf,rank(grotconf));  grotkeep=nets_svds(grotkeep,rank(grotkeep)); 
    beta=pinv([grotconf grotkeep])*y(grot,i); beta(abs(beta)<1e-10)=0; yd(grot,i)=y(grot,i)-grotconf*beta(1:size(grotconf,2));

    if DEMEAN
      yd(grot,i) = nets_demean( yd(grot,i) );
    end
  end

end

