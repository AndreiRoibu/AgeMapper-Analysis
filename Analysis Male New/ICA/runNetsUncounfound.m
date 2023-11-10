addpath FSLNets2

disp('Loading confounds1 ...');
conf12 = h5read('../male_test_conf12.h5','/conf12')';
subjects_conf12 = h5read('../male_test_conf12.h5','/subjects');

disp('Loading brain deltas ...');
load('M_deltas.mat');

fileID = fopen('../../additional_codes/dataset_generation/male_test.txt','r');
formatSpec = '%f';
subjects = fscanf(fileID,formatSpec);
fclose(fileID);
clear fileID formatSpec

X = M;

subjects_to_be_ignored = 21269692;

if size(subjects,1) == size(subjects_conf12,1)
    assert(isequal(subjects, subjects_conf12))
else
    idx_elim = find(subjects == subjects_to_be_ignored);
    subjects(idx_elim, :) = [];
    assert(isequal(subjects, subjects_conf12))
    X(idx_elim, :) = [];
end

% X_deconf = nets_unconfound2(X, conf12);
X_deconf = nets_unconfound(X, conf12);

% Question: should I normalise first and then deconfound? 
% X_deconf and X look quite different

save('M_deltas_deconf.mat','X_deconf','-v7.3')

figure; [counts1, binCenters1] = hist(X(:,1), 100);
[counts2, binCenters2] = hist(X_deconf(:,1), 100);
plot(binCenters1, counts1, 'r-');
hold on;
plot(binCenters2, counts2, 'g-');
grid on;
% Put up legend.
legend1 = sprintf('X');
legend2 = sprintf('X deconf');
legend({legend1, legend2});