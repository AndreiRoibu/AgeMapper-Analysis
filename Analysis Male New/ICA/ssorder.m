function [Aout,Sout] = ssorder(Aref,Ain,Sin);   % reorder columns of Ain and rows of Sin to match Aref
             % it's ok to not include  Sin

N=size(Aref,2);
Nin=size(Ain,2);  % may be less than N

%grot=corrcoef([Aref Ain]);  grot=abs(grot(1:N,N+1:end));
grot=abs(Aref' * Ain);     % using this instead of corrcoef, because the overall mean may be relevant

for ss=1:min(N,Nin)
  grotm=max(grot(:)); [grotii(ss),grotjj(ss)]=find(grot==grotm);
  grot(grotii(ss),:)=0; grot(:,grotjj(ss))=0;
end
[yy,ii]=sort(grotii);
grotjj=grotjj(ii);
Aout=Ain(:,grotjj);

if (exist('Sin'))
  Sout=Sin(grotjj',:);
end

% now fix sign of components
%grot=corrcoef([Aref(:,yy) Aout]); grot=grot(1:Nin,Nin+1:end);
grot=Aref(:,yy)' * Aout;
grot=sign(diag(grot))';

Aout=Aout.*repmat(grot,size(Ain,1),1);
if (exist('Sin'))
  Sout=Sout.*repmat(grot',1,size(Sin,2));
else
  Sout=NaN;
end


