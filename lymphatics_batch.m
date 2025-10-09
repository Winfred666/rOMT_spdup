clear all;
clc;
addpath(genpath('./nii_preprocess'))

% --- CONFIGURATION ---
% Load parameters from the config file, after all preprocess, load psnr_... in set_config_CAA.m
lymph = struct();
lymph.run_Ns = [3,5,7]; % processing steps to run
% Define directories for baseline and source files
lymph.bas_loc_dir = '/data/xym/DEX_MRI/ISO/ISO_52/DCE_nii_baseline';
lymph.src_loc_dir = '/data/xym/DEX_MRI/ISO/ISO_52/DCE_nii_data';

% dst is only for some log files, WARNING: should select an empty folder!
lymph.dst = './temp_output';
% mask for head(bigger) or brain(smaller), used in normalization range + percentage.
% normaliz / percentage will scale all voxels based on the mask region.
lymph.msk = '/data/xym/DEX_MRI/ISO/Template_C57Bl6_n30_brain_ISO_52.nii';

lymph.print = 'n'; % visualize and print image.
lymph.smooth = 0.125 * 2; % smoothing kernel size in mm unit, 2-3 times the voxel size is recommended
lymph.SL = 1; % slice number if want to print image after normalize + smooth


% Helper function to get and unzip files from a directory
function file_list = get_files_from_dir(dir_path)
    if ~exist(dir_path, 'dir')
        error('Directory does not exist: %s', dir_path);
    end
    % delete all the file in folder (do not try for .. and .) except .nii.gz file
    delete_files = dir(fullfile(dir_path, '*'));
    delete_files = delete_files(~contains({delete_files.name}, '.nii.gz'));
    delete_files = delete_files(~ismember({delete_files.name}, {'.', '..'}));
    for i = 1:length(delete_files)
        delete(fullfile(dir_path, delete_files(i).name));
    end
    
    gz_files = dir(fullfile(dir_path, '*.nii.gz'));
    
    file_list = cell(length(gz_files), 1);
    for i = 1:length(gz_files)
        gz_file = fullfile(dir_path, gz_files(i).name);
        fprintf('Unzipping: %s\n', gz_file);
        try
            gunzip(gz_file);
            file_list{i} = erase(gz_file, '.gz'); % Update path in list
        catch ME
            warning('Could not unzip file: %s. Error: %s', gz_file, ME.message);
            file_list{i} = ''; % Mark for removal
        end
    end
    file_list = file_list(~cellfun('isempty', file_list)); % Remove empty entries
end


lymph.bas_loc_files = get_files_from_dir(lymph.bas_loc_dir);
lymph.src_loc_files = get_files_from_dir(lymph.src_loc_dir);
all_files = [lymph.bas_loc_files; lymph.src_loc_files];


% if dst not exist , just create it
if ~exist(lymph.dst, 'dir')
    mkdir(lymph.dst);
    fprintf('Created destination directory: %s\n', lymph.dst);
end


% --------- PROCESSING STEPS ---------

for i = 1:length(lymph.run_Ns)
    switch lower(lymph.run_Ns(i))
        case {1}
            fprintf('1. Format images for Jean Logan \n' );
            lymphatics_JL_images(all_files, lymph.dst, lymph.resol);
        case {2}
            fprintf('2. Multiply image resolution (*v files*)\n' );
            lymph.src_loc = all_files;
            lymph.dst_loc = textread(lymph.dst_list, '%s');
            lymphatics_resize_image(lymph);
        case {3}
            fprintf('3. (Necessary) Realign images with the mean image (*r files*) \n' ); % a must, correct head motions
            lymph.src_loc = all_files;
            lymphatics_realign_image(lymph);
        case {4}
            fprintf('4. Sum images \n' ); % only needed if want to created mean baseline image (however done in 7)
            lymph.src_loc = all_files;
            lymph.dst_loc = lymph.dst;
            lymphatics_sum_image(lymph);
        case {5}
            fprintf('5. (Necessary) Normalize and smooth images (*sn files*)\n' ); % reducing noise and variability in the data.
            lymph.src_loc = all_files;
            % read the output from realigned images, adding 'r' to every filename in src_loc
            for j = 1:length(lymph.src_loc)
                [fdir, fname] = fileparts(lymph.src_loc{j});
                lymph.src_loc{j} = fullfile(fdir, ['r', fname, '.nii']);
            end
            lymphatics_normalize_smooth_images(lymph);
        case {6}
            fprintf('6. Merge all the images to check normalization quality\n' );
            lymph.src_loc = all_files;
            lymph.dst_loc = lymph.dst;
            lymphatics_merge_image(lymph);
        case {7}
            fprintf('7. (Necessary) Convert images to percent images.\n' ); % core step for DCE-MRI analysis, converting the signal intensity changes to percentage changes relative to the baseline
            
            % Process baseline files
            lymph.bas_loc = lymph.bas_loc_files;
            for j = 1:length(lymph.bas_loc)
                [fdir, fname] = fileparts(lymph.bas_loc{j});
                lymph.bas_loc{j} = fullfile(fdir, ['snr', fname, '.nii']);
            end

            % Process source files
            lymph.src_loc = lymph.src_loc_files;
            for j = 1:length(lymph.src_loc)
                [fdir, fname] = fileparts(lymph.src_loc{j});
                lymph.src_loc{j} = fullfile(fdir, ['snr', fname, '.nii']);
            end
            
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


