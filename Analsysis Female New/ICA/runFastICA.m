% This script was adapted from Steve Smith's codes, originally presented in
% https://elifesciences.org/articles/52677, and available here:
% https://www.fmrib.ox.ac.uk/ukbiobank/BrainAgingModes/BrainAgeModesCodes.tar.gz

% This script was adapted by Andrei Roibu, FMRIB, Oxford University,
% March-April 2022, with comments from both the original script and from
% the adapting author.

% We load some required functions
addpath FSLNets FastICA_25 FACS

% Define some useful flags
previously_ran_ICA=0;   % indicate if we have previously ran ICA & saved to disk (1) or not (0); default=0
brain_ages=0;           % indicate if we are using brain ages (1) or deltas (0); default=0
subjects_ica=1;         % indicate if we are loading the subjects-ICA (1) or the features-ICA (0); default=1
only_ICA=0;             % indicate if we're only using ICA (1) or PCA-ICA (0); default=0;
DoPass2=1;              % if previously ran ICA and (1), will repeate DoPass2, else will execute simple code; default=0
load_deconf=1;          % load previously deconfounded values (1) or not (0); default=0;
short_run=1;            % define if we want to run the entire set of J-values (0) or a shorter version (1); default=1


% First we load the subjects x modalities matrix generated with python

if brain_ages==1
    disp('Loading brain ages ...');
    load('M.mat');
else
    disp('Loading brain deltas ...');
    if load_deconf==1
        load('M_deltas_deconf.mat')
        M = X_deconf;
    else
        load('M_deltas.mat')
    end
end

X = M;
clear label

X = nets_normalise(X);

if subjects_ica==1
%     X = nets_normalise(X');
%     X = X';
    C = 2;
    Jvalues = [2:57];
    
else
%     X = nets_normalise(X);
    C = 1;

    if brain_ages==0
        Jvalues = [2:50];
        if short_run==1
            Jvalues = [2:25]
        end
    else
        Jvalues = [2:40];
    end

end

if previously_ran_ICA==0
    disp('DOING ICA...');
    
    % Pre-process the data using Smith et. all implementation, minus the
    % deconfounding and the age-weight scaling
    
%     X=X-nanmedian(X); 
%     grot=nanmedian(abs(X)); 
%     grot(grot<eps)=nanstd(X(:,grot<eps))/1.48; 
%     X=X./grot; X(abs(X)>6)=NaN;
%     X=nets_inormal(X);
%     X(isnan(X))=randn(size(X(isnan(X))))*0.01;
%     X=nets_demean(X);
%     clear grot


    % ICA on IDPs, including dimensionality estimation

    if subjects_ica==1
        disp('Doing FastICA in subject direction ...');
        disp(size(X'));
        if only_ICA==0
            disp('Doing PCA-ICA...');
          
            if brain_ages==0
                output_name_subject_ICA = 'PcaIca_subject_deltas_transpose_steve_norm.mat';
                if load_deconf==1
                    output_name_subject_ICA = 'PcaIca_subject_deltas_transpose_steve_norm_deconf.mat';
                    if short_run==1
                        output_name_subject_ICA = 'PcaIca_subject_deltas_transpose_steve_norm_deconf_short.mat';
                    end
                end
            else
                output_name_subject_ICA = 'PcaIca_subject_transpose_steve_norm_ages.mat';
            end

            [icaS,icaA,pcaU,pcaS,pcaV,correlations_dict]=ICAdim2(X', Jvalues, 2, 10, 3000); 

            save(output_name_subject_ICA,'pca*','ica*','correlations_dict');

        else
            disp('Only doing ICA...');
            [icaS,icaA,correlations_dict]=ICAdim2_noPCA(X', 2); 
            save(output_name_subject_ICA,'ica*','correlations_dict');
        end
        tempdir = 'subjects_ica_figs';
        mkdir subjects_ica_figs
    else
        disp('Doing FastICA in modality direction ...');
        disp(size(X));
        if brain_ages==0
            
            output_name_features_ICA = 'PcaIca_features_deltas_10000_x_100_x_0.9_deltas.mat';
            if load_deconf==1
                output_name_features_ICA = 'PcaIca_features_deltas_10000_x_100_x_0.9_deltas_deconf.mat';
                if short_run==1
                    output_name_features_ICA = 'PcaIca_features_deltas_10000_x_100_x_0.9_deltas_deconf_short.mat';
                end
            end

            [icaS,icaA,pcaU,pcaS,pcaV,correlations_dict]=ICAdim2(X, Jvalues, 1, 100, 10000);

        else
            [icaS,icaA,pcaU,pcaS,pcaV,correlations_dict]=ICAdim2(X, Jvalues, 1, 100, 10000);
            output_name_features_ICA = 'PcaIca_features_deltas_10000_x_100_x_0.9_ages.mat';
        end

        save(output_name_features_ICA,'pca*','ica*','correlations_dict');
        tempdir = 'features_ica_figs';
        mkdir features_ica_figs;
    end

    make_corr_plots(correlations_dict, 0.8);
    make_corr_plots(correlations_dict, 0.7);
    make_corr_plots(correlations_dict, 0.5);

    FolderName = tempdir;   % Your destination folder
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = num2str(get(FigHandle, 'Number'));
      set(0, 'CurrentFigure', FigHandle);
      savefig(fullfile(FolderName, [FigName '.fig']));
    end

    if DoPass2==1

        disp('ICA ALREADY DONE! DOING PASS 2 NOW....');
        if subjects_ica==1
            if brain_ages==0
                load_matrix_name = 'PcaIca_subject_deltas_transpose_steve_norm.mat';
                if load_deconf==1
                    load_matrix_name = 'PcaIca_subject_deltas_transpose_steve_norm_deconf.mat';
                    if short_run==1
                        load_matrix_name = 'PcaIca_subject_deltas_transpose_steve_norm_deconf_short.mat';
                    end
                end
            else
                load_matrix_name = 'PcaIca_subject_transpose_steve_norm_ages.mat';
            end
        else
            if brain_ages==0
                load_matrix_name = 'PcaIca_features_deltas_10000_x_100_x_0.9_deltas.mat';
                if load_deconf==1
                    load_matrix_name = 'PcaIca_features_deltas_10000_x_100_x_0.9_deltas_deconf.mat';
                    if short_run==1
                        load_matrix_name = 'PcaIca_features_deltas_10000_x_100_x_0.9_deltas_deconf_short.mat';
                    end
                end
            else
                load_matrix_name = 'PcaIca_features_deltas_10000_x_100_x_0.9_ages.mat';
            end
        end
    
        disp("Loading: "), disp(load_matrix_name)
        load(load_matrix_name)
    
        if subjects_ica==1
            [icaS,icaA,pcaU,pcaS,pcaV] = ICAdim2_Pass2(X', Jvalues, correlations_dict, C);
        else
            [icaS,icaA,pcaU,pcaS,pcaV] = ICAdim2_Pass2(X, Jvalues, correlations_dict, C);
        end
        output_file = strcat(load_matrix_name(1:size(load_matrix_name,2)-4) , '_pass2.mat');
        save(output_file,'pca*','ica*','correlations_dict');


    end

else

    disp('ICA ALREADY DONE...');
    if subjects_ica==1
        if brain_ages==0
            load_matrix_name = 'PcaIca_subject_deltas_transpose_steve_norm.mat';
            if load_deconf==1
                load_matrix_name = 'PcaIca_subject_deltas_transpose_steve_norm_deconf.mat';
            end
        else
            load_matrix_name = 'PcaIca_subject_transpose_steve_norm_ages.mat';
        end
    else
        if brain_ages==0
            load_matrix_name = 'PcaIca_features_deltas_10000_x_100_x_0.9_deltas.mat';
            if load_deconf==1
                load_matrix_name = 'PcaIca_features_deltas_10000_x_100_x_0.9_deltas_deconf.mat';
            end
        else
            load_matrix_name = 'PcaIca_features_deltas_10000_x_100_x_0.9_ages.mat';
        end
    end

    disp("Loading: "), disp(load_matrix_name)
    load(load_matrix_name)

    if DoPass2==1
        if subjects_ica==1
            [icaS,icaA,pcaU,pcaS,pcaV] = ICAdim2_Pass2(X', Jvalues, correlations_dict, C);
        else
            [icaS,icaA,pcaU,pcaS,pcaV] = ICAdim2_Pass2(X, Jvalues, correlations_dict, C);
        end
        output_file = strcat(load_matrix_name(1:size(load_matrix_name,2)-4) , '_pass2.mat');
        save(output_file,'pca*','ica*','correlations_dict');
    else
        vs = values(correlations_dict); ks=keys(correlations_dict);
        vs1 = []; ks1 = [];
        idx=1; CorrThresh=0.8;
        for v=vs
        ks1 = [ks1  str2double(cell2mat(ks(idx)))];
        idx=idx+1;
        vs1 = [vs1 sum(mean(cell2mat(v))>CorrThresh)];
        end
        [ks1_sorted, a_order] = sort(ks1);
        vs1_sorted = vs1(a_order);
        figure; plot(ks1_sorted,vs1_sorted); title(sprintf('Number of Strong correlations at thr=%d',CorrThresh)); drawnow;
    end
end