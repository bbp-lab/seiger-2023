function []=AN_RSA_Searchlight(subj,mask_name)

% Start Date: 3/31/22
% Contact: Nicco Reggente, Ph.D. (nicco@advancedconsciousness.org)

tic;
disp(['Running Subject....' subj])
%% Hard Code
flags.default_size=[91 109 91];
%Of note is that there should be 120 volumes in total; the first 60 are
%anxiety words and the last 60 are neutral words.
%We will hardcode Anxiety as 1 and Neutral as 2;
trial_labels=[ones(1,60),ones(1,60)*2];
[trial_names{1:60}]=deal('Anxiety');
[trial_names{61:120}]=deal('Neutral');
regressors={'motion'};

%% Leverage Naming Conventions
naming.pre_subj='sub_'; %The standard prefix on files. Include underscores if they are used.
naming.post_subj='_proc'; %The standard suffix on files. Include underscores if they are used.
naming.subj_folder='sub-';

%% Define Paths
paths.top='D:\PROJECTS\AN_RSA\';
paths.reference=[paths.top 'Reference/'];
paths.data=[paths.top 'Data/' naming.subj_folder subj '/'];
paths.masks=[paths.top 'Masks/'];
paths.ml='D:\PROJECTS\MATLAB_PATH'; addpath(genpath(paths.ml)); %add toolboxes to the path
paths.save=[paths.top 'Subject_Results/Searchlight/'];

%% Get Patterns of Interest
%Normally...use the script 'create_patterns' to get patterns for each data type for
%each mask. Since we are only getting 2 trial types, we can do it all in one function call.

%First step is to index the total available trials with "active" trials. In
%this particular use case, we are using excess motion to determine

active_trials=zeros(1,size(trial_labels,2));
for i=1:size(regressors,1)
    temp_active_trials=dlmread([paths.data naming.pre_subj subj '_' regressors{i} '.txt']);
    active_trials=active_trials+temp_active_trials;
end

active_trials_idx=find(active_trials==0); %0 means they are not excluded.
num_trials_excluded=length(find(active_trials==1)); %Potentially useful down the road. Will be written to results.

trial_types=unique(trial_labels);

patterns=cell(1,length(trial_types));

for t=1:length(trial_types)
    trials_of_interest=intersect(find(trial_labels==trial_types(t)),active_trials_idx); %Find the indices that for trials of interest AND active.
    [patterns{t}]=create_SL_patterns(subj,paths,mask_name, trials_of_interest, flags, naming);
end
clear temp_patterns


%% Create Across Maps
disp('Creating Across Patterns...')
for p=1:(size(patterns,2)-1);%Loop over each type of pair.
    Temp1=patterns{p};
    Temp2=patterns{p+1};
    %     ppm=ParforProgressbar(size(Temp1,2),'showWorkerProgress',true,'progressBarUpdatePeriod',4,'title','Across Pattern Building');
    for sl=1:size(patterns{p},2)% Loop over searchlights
        A=size(Temp1{sl},1);
        B=size(Temp2{sl},1);
        C=zeros(1,A*B);
        for a=1:A
            for b=1:B
                C(1,(a*b))=corr(cell2mat(Temp1{sl}(a,:))',cell2mat(Temp2{sl}(b,:))');
            end
        end
        temp_r{sl}=mean(atanh(C));%Grab the non-zero values and do an r-to-z transform.
        %         ppm.increment()
        progress(sl,size(patterns{p},2),20);
    end
    %     delete(ppm);
    Across_R_Values{p}=cell2mat(temp_r);
    clear temp_r
end


%% Create Within Maps
disp('Creating Within Patterns...')

for p=1:size(patterns,2)%Loop over each type of pattern.
    %     ppm=ParforProgressbar(size(patterns{p},2),'showWorkerProgress',true,'progressBarUpdatePeriod',4,'title','Within Pattern Building');
    for sl=1:size(patterns{p},2)% Loop over searchlights
        A=tril(corr(cell2mat(patterns{p}{sl})'),-1);%Within each searchlight sphere, correlate all the exmemplars
        temp_r{sl}=mean(atanh(A(find(A))));%Grab the non-zero values and do an r-to-z transform.
        %         ppm.increment();
        progress(sl,size(patterns{p},2),20);
    end
    Within_R_Values{p}=cell2mat(temp_r);
    %     delete(ppm);
    clear temp_r

end


%% Create Within vs. Across Comparisons
%Statistically needs to be done at the group level, but can create the
%subtraction maps here.
across_map=Across_R_Values{1};
within_anxiety_map=Within_R_Values{1};
within_neutral_map=Within_R_Values{2};

within_anxiety_map_minus_within_neutral_map=within_anxiety_map-within_neutral_map;
within_neutral_map_minus_within_anxiety_map=within_neutral_map-within_anxiety_map;

across_minus_within_anxiety_map=across_map-within_anxiety_map;
across_minus_within_neutral_map=across_map-within_neutral_map;

within_anxiety_minus_across=within_anxiety_map-across_map;
within_neutral_minus_across=within_neutral_map-across_map;


%% Write Outputs

V=spm_vol([paths.masks '/' mask_name '.nii']);
temp_vols=spm_read_vols(V);
active_mask_idx=find(temp_vols);
V.dt=[16 0];

%Write Within Me Maps
vols=temp_vols;
vols(active_mask_idx)=within_anxiety_map; %assign R values to the active voxels.
V.fname=[paths.save subj '_' mask_name '_Searchlight_Within_Anxiety.nii'];
spm_write_vol(V,vols);

%Write Within NotMe Maps
vols=temp_vols;
vols(active_mask_idx)=within_neutral_map; %assign R values to the active voxels.
V.fname=[paths.save subj '_' mask_name '_Searchlight_Within_Neutral.nii'];
spm_write_vol(V,vols);

%Write Across Maps
vols=temp_vols;
vols(active_mask_idx)=across_map; %assign R values to the active voxels.
V.fname=[paths.save subj '_' mask_name '_Searchlight_Across.nii'];
spm_write_vol(V,vols);

%Write Subtraction Maps
vols=temp_vols;
vols(active_mask_idx)=within_anxiety_map_minus_within_neutral_map; %assign R values to the active voxels.
V.fname=[paths.save subj '_' mask_name '_Searchlight_Within_Anxiety_Minus_Within_Neutral.nii'];
spm_write_vol(V,vols);

vols=temp_vols;
vols(active_mask_idx)=within_neutral_map_minus_within_anxiety_map; %assign R values to the active voxels.
V.fname=[paths.save subj '_' mask_name '_Searchlight_Within_Neutral_Minus_Within_Anxiety.nii'];
spm_write_vol(V,vols);

vols=temp_vols;
vols(active_mask_idx)=across_minus_within_anxiety_map; %assign R values to the active voxels.
V.fname=[paths.save subj '_' mask_name '_Searchlight_Across_Minus_Within_Anxiety.nii'];
spm_write_vol(V,vols);

vols=temp_vols;
vols(active_mask_idx)=across_minus_within_neutral_map; %assign R values to the active voxels.
V.fname=[paths.save subj '_' mask_name '_Searchlight_Across_Minus_Within_Neutral.nii'];
spm_write_vol(V,vols);

vols=temp_vols;
vols(active_mask_idx)=within_anxiety_minus_across; %assign R values to the active voxels.
V.fname=[paths.save subj '_' mask_name '_Searchlight_Within_Anxiety_Minus_Across.nii'];
spm_write_vol(V,vols);

vols=temp_vols;
vols(active_mask_idx)=within_neutral_minus_across; %assign R values to the active voxels.
V.fname=[paths.save subj '_' mask_name '_Searchlight_Within_Neutral_Minus_Across.nii'];
spm_write_vol(V,vols);

disp(['Done with Subject ' subj '. That Took: ' num2str(toc/60) ' minutes.'])

end