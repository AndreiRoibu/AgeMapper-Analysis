%
% nets_hierarchy - create hierarchical clustering figure
% Steve Smith, 2012-2014
%
% [hier_order,linkages] = nets_hierarchy(netmat)
% [hier_order,linkages] = nets_hierarchy(netmatL,netmatH)
% [hier_order,linkages] = nets_hierarchy(netmatL,netmatH,DD,sumpics)
% [hier_order,linkages] = nets_hierarchy(netmatL,netmatH,DD,sumpics,colour_threshold)
%
%   netmatL is the net matrix shown below the diagonal, and drives the hierarchical clustering.
%           It is typically the z-stat from the one-group-t-test across all subjects' netmats
%   netmatH is the net matrix shown above the diagonal, for example partial correlation.
%   DD is the list of good nodes (needed to map the current set of nodes back to originals), e.g., ts.DD
%   sumpics is the name of the directory containing summary pictures for each component, without the .sum suffix
%   hier_order (output) is the index numbers of the good nodes, as reordered by the hierarchical clustering.
%   colour_threshold (default 0.75) specifies at what level the dendrogram colouring stops
%

function [dpRSN,yyRSN] = nets_hierarchy_andrei_mod(netmatL,varargin);    %%%% hierarchical clustering figure

netmatH=netmatL;
if nargin>1
  netmatH=varargin{1};
end

DD=[];
if nargin>2
  DD=varargin{2};
end

sumpics='';
if nargin>3
  sumpics=varargin{3};
end

colour_threshold=0.75;
if nargin>4
  colour_threshold=varargin{4};
end

% labels = ["T1 nonlinear","T1 linear","jacobian","vbm","T2 nonlinear","T2 lesions","swi","rsfmri 0","rsfmri 1","rsfmri 2","rsfmri 3","rsfmri 4","rsfmri 5","rsfmri 6","rsfmri 7","rsfmri 8","rsfmri 9","rsfmri 10","rsfmri 11","rsfmri 12","rsfmri 13","rsfmri 14","rsfmri 15","rsfmri 16","rsfmri 17","rsfmri 18","rsfmri 19","rsfmri 20","rsfmri 21","rsfmri 22","rsfmri 23","rsfmri 24","tfmri 1","tfmri 2","tfmri 5","tfmri c 1","tfmri c 2","tfmri c 5","tracts","tbss FA s","tbss ICVF s","tbss ISOVF s","tbss L1 s","tbss L2 s","tbss L3 s","tbss MD s","tbss MO s","tbss OD s","tbss FA","tbss ICVF","tbss ISOVF","tbss L1","tbss L2","tbss L3","tbss MD","tbss MO","tbss OD"];
labels = ["T1 Nonlinear", "T1 Linear", "Jacobian", "VBM", "T2 Nonlinear", "T2 Lessions", "SWI", "rsfMRI-0", "rsfMRI-1","rsfMRI-2", "rsfMRI-3", "rsfMRI-4", "rsfMRI-5", "rsfMRI-6", "rsfMRI-7", "rsfMRI-8", "rsfMRI-9", "rsfMRI-10","rsfMRI-11", "rsfMRI-12", "rsfMRI-13", "rsfMRI-14", "rsfMRI-15", "rsfMRI-16", "rsfMRI-17", "rsfMRI-18", "rsfMRI-19","rsfMRI-20", "rsfMRI-21", "rsfMRI-22", "rsfMRI-23", "rsfMRI-24", "tfMRI-1", "tfMRI-2", "tfMRI-5", "tfMRI-COPE-1","tfMRI-COPE-2", "tfMRI-COPE-5", "Summed Tracts", "TBSS FA", "TBSS ICVF", "TBSS ISOVF", "TBSS L1", "TBSS L2","TBSS L3", "TBSS MD", "TBSS MO", "TBSS OD", "FA", "ICVF", "ISOVF", "L1", "L2", "L3", "MD", "MO", "OD"]; 
 
%grot=prctile(abs(netmatL(:)),99); sprintf('99th%% abs value below diagonal is %f',grot)
%netmatL=netmatL/grot;
%grot=prctile(abs(netmatH(:)),99); sprintf('99th%% abs value above diagonal is %f',grot)
%netmatH=netmatH/grot;

grot1=prctile(abs(netmatL(triu(ones(size(netmatL)),1)>0)),99); sprintf('99th%% abs value below diagonal is %f',grot1)
grot2=prctile(abs(netmatH(triu(ones(size(netmatH)),1)>0)),99); sprintf('99th%% abs value above diagonal is %f',grot2)
grot=grot1;
netmatL=netmatL/grot;
netmatH=netmatH/grot;

usenet=netmatL;
usenet(usenet<0)=0;   % zero negative entries.....seems to give nicer hierarchies

H1=0.76; H2=0.8;    % sub-image sizes for case of no summary pics
[grotA,grotB,grotC]=fileparts(sumpics); if size(grotA,1)==0, grotA='.'; end; sumpics=sprintf('%s/%s.sum',grotA,grotB);
if exist(sumpics) > 0
  H1=0.73; H2=0.8;    % sub-image sizes for default of 1-slice summary pics
  grot=imread(sprintf('%s/0000.png',sumpics));  % read the first summary image to find out how many slices are in each
  if abs((size(grot,1)/size(grot,2)-1))>0.5
    H1=0.6; H2=0.8;   % sub-image sizes for 3-slice summary pics
  end
end

figure('position',[10 10 1000 500]);  clear y;  N=size(netmatL,1);  gap=.5/(N+1);  
grot=prctile(abs(usenet(:)),99); usenet=max(min(usenet/grot,1),-1)/2;
for J = 1:N, for I = 1:J-1,   y((I-1)*(N-I/2)+J-I) = 0.5 - usenet(I,J);  end; end;
  yyRSN=linkage(y,'ward');
subplot('position',[0 H2 1 1-H2]);
  if exist('octave_config_info')~=0
    subplot('position',[gap H2 1-2*gap 1-H2]);
  end;
  [dh,dt,dpRSN]=dendrogram(yyRSN,0,'colorthreshold',colour_threshold,'Labels',labels); set(gca,'ytick',[],'TickLength',[0,0]);
  if exist('octave_config_info')==0
    set(dh,'LineWidth',3);
  end;
subplot('position',[gap 0 1-2*gap H1-0.01]); i=dpRSN;
  grot=tril(netmatL(i,i)) + triu(netmatH(i,i)); grot=max(min(grot,0.95),-0.95);  grot(eye(length(grot))>0)=Inf;
  colormap('jet');
  grotc=colormap;  grotc(end,:)=[.8 .8 .8];  colormap(grotc);  imagesc(grot,[-1 1]); axis off; daspect('auto');
grott=sprintf('%s.png',tempname);
  if length(DD)>0
    system(sprintf('slices_summary %s %s %s',sumpics,grott,num2str(DD(dpRSN)-1)));
  end
  if exist(grott) == 2
    grot=imread(grott);  subplot('position',[gap H1 1-2*gap H2-H1-0.035]); imagesc(grot); axis off; system(sprintf('/bin/rm %s*',grott)); daspect('auto');
  end;
set(gcf,'PaperPositionMode','auto');

%print('-dpng',sprintf('hier.png'));

