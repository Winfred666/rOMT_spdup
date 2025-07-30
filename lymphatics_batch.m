clear all;
clc;
addpath(genpath('./nii_preprocess'))

% --- CONFIGURATION ---
% Load parameters from the config file, after all preprocess, load psnr_... in set_config_CAA.m
lymph = struct();
lymph.run_Ns = [3,5,7]; % processing steps to run
% lymph.src = '/home/xym/Desktop/MRI_ROMT/rOMT_spdup-main/data/ours/test1/total_src.txt'; % source image
dataset = 'DEXI';
dataset_num = '084';
lymph.src = sprintf('/data/xym/DEX_MRI/%s/%s_%s/src_total.txt', dataset,dataset,dataset_num);
% dst is only for some log files.
% lymph.dst = '/home/xym/Desktop/MRI_ROMT/rOMT_spdup-main/data/ours/test1/DCE_nii_preprocess_log'; % destination directory
lymph.dst = sprintf('/data/xym/DEX_MRI/%s/%s_%s/DCE_nii_preprocess_log', dataset, dataset, dataset_num);
% mask for brain (still the same), used in normalization range.
lymph.msk = sprintf('/data/xym/DEX_MRI/%s/Template_C57Bl6_n30_brain_%s_%s.nii', dataset, dataset, dataset_num);

lymph.print = 'n'; % visualize and print image.
lymph.smooth = 0.125 * 4; % smoothing kernel size in mm unit, 2-3 times the voxel size is recommended
lymph.SL = 1; % slice number if want to print image after normalize + smooth







% if there is no *.nii but only *.nii.gz files in src, use gunzip to unpack all *.gz
fileID = fopen(lymph.src, 'r');
if fileID == -1
    error('Cannot open source file: %s', lymph.src);
end
file_list = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);
file_list = file_list{1};

file_updated = false; % Flag to check if the file needs to be rewritten

% Iterate over the list of files to find and unzip .nii.gz
for k = 1:length(file_list)
    current_file = file_list{k};
    if endsWith(current_file, '.nii.gz')
        fprintf('Unzipping: %s\n', current_file);
        try
            gunzip(current_file);
            file_list{k} = erase(current_file, '.gz'); % Update path in list
            file_updated = true;
        catch ME
            warning('Could not unzip file: %s. Error: %s', current_file, ME.message);
        end
    end
end

% if dst not exist , just create it
if ~exist(lymph.dst, 'dir')
    mkdir(lymph.dst);
    fprintf('Created destination directory: %s\n', lymph.dst);
end

% If any files were unzipped, rewrite the source file with the updated .nii paths
if file_updated
    fprintf('Rewriting source file with updated .nii paths: %s\n', lymph.src);
    fileID = fopen(lymph.src, 'w');
    if fileID == -1
        error('Cannot open source file for writing: %s', lymph.src);
    end
    fprintf(fileID, '%s\n', file_list{:});
    fclose(fileID);
    fprintf('Source file updated successfully.\n');
end


% --------- PROCESSING STEPS ---------

for i = 1:length(lymph.run_Ns)
    switch lower(lymph.run_Ns(i))
        case {1}
            fprintf('1. Format images for Jean Logan \n' );
            lymphatics_JL_images(lymph.src, lymph.dst, lymph.resol);
        case {2}
            fprintf('2. Multiply image resolution (*v files*)\n' );
            lymph.src_loc = textread(lymph.src, '%s');
            lymph.dst_loc = textread(lymph.dst_list, '%s');
            lymphatics_resize_image(lymph);
        case {3}
            fprintf('3. (Necessary) Realign images with the mean image (*r files*) \n' ); % a must, correct head motions
            lymph.src_loc = textread(lymph.src, '%s');
            lymphatics_realign_image(lymph);
        case {4}
            fprintf('4. Sum images \n' ); % only needed if want to created mean baseline image (however done in 7)
            lymph.src_loc = textread(lymph.src, '%s');
            lymph.dst_loc = lymph.dst;
            lymphatics_sum_image(lymph);
        case {5}
            fprintf('5. (Necessary) Normalize and smooth images (*sn files*)\n' ); % reducing noise and variability in the data.
            lymph.src_loc = textread(lymph.src, '%s');
            % read the output from realigned images, adding 'r' to every filename in src_loc
            for j = 1:length(lymph.src_loc)
                [fdir, fname] = fileparts(lymph.src_loc{j});
                lymph.src_loc{j} = fullfile(fdir, ['r', fname, '.nii']);
            end
            lymphatics_normalize_smooth_images(lymph);
        case {6}
            fprintf('6. Merge all the images to check normalization quality\n' );
            lymph.src_loc = textread(lymph.src, '%s');
            lymph.dst_loc = lymph.dst;
            lymphatics_merge_image(lymph);
        case {7}
            fprintf('7. (Necessary) Convert images to percent images.\n' ); % core step for DCE-MRI analysis, converting the signal intensity changes to percentage changes relative to the baseline
            lymph.src_loc = textread(lymph.src, '%s');

            % read output from smoothed images, add 'sn' to every filename in src_loc
            for j = 1:length(lymph.src_loc)
                [fdir, fname] = fileparts(lymph.src_loc{j});
                lymph.src_loc{j} = fullfile(fdir, ['snr', fname, '.nii']);
            end

            % filter files in src_loc list that include "baseline" into bas_loc
            baseline_idxs = contains(lymph.src_loc, 'baseline');
            
            lymph.bas_loc = lymph.src_loc(baseline_idxs);
            lymph.src_loc = lymph.src_loc(~baseline_idxs);
            
            % lymph.bas_loc = textread(lymph.bas, '%s');
            lymphatics_percent_image(lymph);
        case {8}
            fprintf('8. Re-format images for saggital orientation.\n' );
            if lymph.modal == 1
                lymph.ref = '/usr/local/matcodes/spm12_batch/Lymphatics_Prj/070912_template/pbase_snrv_Gadospin_062912A_E40.img';
            elseif lymph.modal == 2
                lymph.ref = '/usr/local/matcodes/spm12_batch/Lymphatics_Prj/032613_baby_template/rrpbase_snrv_Magnevist_23_E65_032613.img';
            elseif lymph.modal == 3
                lymph.ref = '/usr/local/matcodes/spm12_batch/Lymphatics_Prj/071714_Prone_template/rpbase_snrv_Magnevist_40_water_Cistern_prone_E50_071714A.img';
            elseif lymph.modal == 4
                lymph.ref = '/usr/local/matcodes/spm12_batch/Lymphatics_Prj/082014_Lateral_template/rpbase_snrv_Magnevist_40_water_Cistern_lateral_E71_082014A.img';
            elseif lymph.modal == 5
                lymph.ref = '/usr/local/matcodes/spm12_batch/Lymphatics_Prj/072815_supine_template/rpbase_snrv_Magnevist_40_water_cistern_supine_NIDA_E55_072115B.nii';
            end
            
            [lymph.dst, ~] = fileparts(lymph.src);
            lymph.oth_loc = textread(lymph.oth, '%s');
            lymphatics_realign_sagittal_image(lymph);
        case {9}  
            fprintf('9. Normalize brain.\n' ); % usually warp ATLS, this individual warp is for certain brain damaged rats or mouse.
            lymphatics_warp_batch;
        otherwise
            disp('Unknown option.')
    end
end


