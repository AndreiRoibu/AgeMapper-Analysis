
function [LoadedMatched] = nets_load_match(infile,to_match,varargin);

DELETE_COLUMNS=1;
if nargin==3
  DELETE_COLUMNS=0;
end

grot=load(infile);

if length(grot(:,1)) > length(unique(grot(:,1)))
  disp('Warning - looks like duplicate subject IDs in input file. Processing anyway but may be wrong.');
end

LoadedMatched = zeros(size(to_match,1),size(grot,2)-DELETE_COLUMNS) / 0;
[~,grotI,grotJ]=intersect(grot(:,1),to_match(:,1));
LoadedMatched(grotJ,:)=grot(grotI,(1+DELETE_COLUMNS):end);

