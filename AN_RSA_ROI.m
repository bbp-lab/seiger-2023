function []=AN_RSA_ROI(subj,mask_name)

% Start Date: 3/30/22
% Contact: Nicco Reggente, Ph.D. (nicco@advancedconsciousness.org)

%Workspace Running
%AN_RSA_ROI('4011','Body_Localizer')


%% Input Insights
%Inputs must take the following form and [var type]:

%data_types [cell array of numbers] = Which data types (max 2) will be used
%to run the RSA comparison of within vs. between. (e.g. {1,2}

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
paths.reference=[paths.top 'Reference'];
paths.data=[paths.top 'Data/' naming.subj_folder subj '/'];
paths.masks=[paths.top 'Masks/'];
paths.ml='D:\PROJECTS\MATLAB_PATH'; addpath(genpath(paths.ml)); %add toolboxes to the path
paths.save=[paths.top 'Subject_Results/ROIs/'];

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

for t=1:length(trial_types)
    trials_of_interest=intersect(find(trial_labels==trial_types(t)),active_trials_idx); %Find the indices that for trials of interest AND active.
    [patterns{t}]=create_ROI_patterns(subj,paths,mask_name, trials_of_interest, flags, naming);
end

%% Create Across Values
disp('Creating Across Patterns...')
for p=1:(size(patterns,2)-1)%Loop over each type of pair.
    A=size(patterns{p},1);
    B=size(patterns{p+1},1);
    C=zeros(1,A*B);
    for a=1:A
        for b=1:B
            C(1,(a*b))=corr(cell2mat(patterns{p}(a,:))',cell2mat(patterns{p+1}(b,:))');
        end
        progress(a,A,20);
    end
    temp_r=mean(atanh(C(find(C))));%Grab the non-zero values and do an r-to-z transform.
    Across_R_Values{p}=temp_r;
    clear temp_r
end

%% Create Within Values
disp('Creating Within Patterns...')
for p=1:size(patterns,2)%Loop over each type of pattern.
    A=tril(corr(cell2mat(patterns{p})'),-1);%Within ROI, correlate all the exmemplars
    temp_r=mean(atanh(A(find(A))));%Grab the non-zero values and do an r-to-z transform.
    Within_R_Values{p}=temp_r;
end

clear temp_r


%% Write Outputs

header={'Within_Anxiety','Within_Neutral','Across','Trials Excluded'};
data=[Within_R_Values{1},Within_R_Values{2},Across_R_Values{1}, num_trials_excluded];
save_file=[paths.save subj '_' mask_name '.txt'];

save_data_with_headers(header,data,save_file);


end