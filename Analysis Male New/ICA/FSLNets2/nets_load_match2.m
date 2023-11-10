
function [LoadedMatched] = nets_load_match(infile,to_match);

grot=load(infile);
LoadedMatched = zeros(size(to_match,1),size(grot,2)-2) / 0;
grotID=grot(:,2)*1e7+grot(:,1);
[~,grotI,grotJ]=intersect(grotID,to_match(:,1));
LoadedMatched(grotJ,:)=grot(grotI,3:end);

