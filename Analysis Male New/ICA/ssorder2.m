function [AoutRef,SoutRef,Aout1,Sout1,Aout2,Sout2] = ssorder2(AinRef,SinRef,Ain1,Sin1,Ain2,Sin2); % re-order in1 and in2 to match inRef, then re-order all 3 according to in1 vs in2 matching quality

[Aout1,Sout1]=ssorder(AinRef,Ain1,Sin1);
[Aout2,Sout2]=ssorder(AinRef,Ain2,Sin2);

ii=size(Aout1,2);
grot=corrcoef([Aout1 Aout2]);
grot=diag(grot(1:ii,ii+1:end));

[yy,ii]=sort(grot,'descend');

AoutRef=AinRef(:,ii);
SoutRef=SinRef(ii,:);
Aout1=Aout1(:,ii);
Sout1=Sout1(ii,:);
Aout2=Aout2(:,ii);
Sout2=Sout2(ii,:);

