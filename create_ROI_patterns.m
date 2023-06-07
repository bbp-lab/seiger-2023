function [patterns_of_interest]=create_ROI_patterns(subj,paths, mask_name, trials_of_interest, flags, naming)

% Start Date: 3/30/22
% Contact: Nicco Reggente, Ph.D. (nicco@advancedconsciousness.org)

%This script will output a matrix with rows of features (i.e. voxel values)
%with columns of instances (i.e. trials / exemplars).
% This script is specific to an ROI procedure

%% Input Insights
%Inputs must take the following form and [var type]:

%subj [string] = Subject Number
%mask [vector] = Indices for active voxels (must be in same space as image)
%paths [struct] = Paths structure with relevant directory paths
%naming [struct] = Naming conventions for loading in subject data
%trials_of_interest [vector] = Indices for active trials to be loaded in
%from a subject's 4D file.

%% Bring in Subject-Specific Data

%Load in Subject's 4D file (containing all trial types)
subj_4D_file=[paths.data naming.pre_subj subj naming.post_subj '.nii'];
if ~exist(subj_4D_file)
    gunzip([subj_4D_file '.gz']); %Maybe was not gunzipped first.
end
V=spm_vol(subj_4D_file);
vols=spm_read_vols(V);
%only bring in beta maps that represent trials of interest
vols_of_interest=vols(:,:,:,trials_of_interest);
clear vols V %Keep tabs on memory to run smoothly

%% Mask the Data with voxels in a mask
V=spm_vol([paths.masks mask_name '.nii']);
mask_vols=spm_read_vols(V);
mask_indices=find(mask_vols);

disp('Creating Patterns for ROI...')
for t=1:length(trials_of_interest)
    for m=1:length(mask_indices)
        [a, b, c]=ind2sub(flags.default_size, mask_indices(m));
        patterns_of_interest{t,m}=reshape_to_1(vols_of_interest(a,b,c,t));
    end
    progress(t,length(trials_of_interest),10);
end

end