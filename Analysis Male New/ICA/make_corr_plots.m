function [] = make_corr_plots(correlations_dict, CorrThresh)
%MAKE_CORR_PLOTS Summary of this function goes here
%   Detailed explanation goes here

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

end

