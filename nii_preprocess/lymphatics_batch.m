clear all;
clc;

fprintf('1. Format images for Jean Logan \n' );
fprintf('2. Multiply image resolution (*v files*)\n' );
fprintf('3. Realign images with the mean image (*r files*) \n' );
fprintf('4. Sum images \n' );
fprintf('5. Normalize and smooth images (*sn files*)\n' );
fprintf('6. Merge all the images to check normalization quality\n' );
fprintf('7. Convert images to percent images.\n' );
fprintf('8. Re-format images for saggital orientation.\n' );
fprintf('9. Normlize brain.\n' );

run_N=input('Enter the numer:');

switch lower(run_N)
    case {1}
        lymph.src=input('Enter the source directory: ','s');
        lymph.dst=input('Enter the destination directory: ','s');
        lymph.resol=input('Cut image resolution by half ? (ex:y or n): ','s');
        lymphatics_JL_images(lymph.src,lymph.dst,lymph.resol);
    case {2}
        lymph.species=input('Enter the species 1.Rat (x10) 2.Mouse (x20): '); %4/3/18
        lymph.src=input('Enter the source files(ex:/060112A/Magnevist_060112A_E9.nii): ','s');
        lymph.dst=input('Enter the output files(ex:/060112A/v_Magnevist_060112A_E9.nii): ','s');
        lymph.src_loc=textread(lymph.src,'%s');
        lymph.dst_loc=textread(lymph.dst,'%s');
        lymphatics_resize_image(lymph);
    case {3}
        lymph.src=input('Enter the src file (always align all image with the mean): ','s');
        lymph.dst=input('Enter the destination directory to save the batch file: ','s');
        lymph.src_loc=textread(lymph.src,'%s');
        lymphatics_realign_image(lymph);
    case {4}
        lymph.src=input('Enter the source files(ex:src.txt): ','s');
        lymph.dst=input('Enter the output file name(ex:/060112A/mean.nii): ','s');
        lymph.src_loc=textread(lymph.src,'%s');
        lymph.dst_loc=lymph.dst;
        lymphatics_sum_image(lymph);
    case {5}
        lymph.src=input('Enter the source files(ex:list.txt): ','s');
        lymph.msk=input('Enter the mask image for NORMALIZATION(ex:img_msk.nii): ','s');
        lymph.smooth=input('Enter smoothing kernel in mm (ex:1): ');
        lymph.dst=input('Enter the destination directory(ex:/060112A): ','s');
        %Added on 03/20/17
        lymph.print=input('Save output as png ?(y or n): ','s');
        if strcmp(lymph.print,'y')
            lymph.SL=input('Enter the slice number for print out(135 works well): ');
        end
        lymph.src_loc=textread(lymph.src,'%s');
        lymphatics_normalize_smooth_images(lymph);
    case {6}
        lymph.src=input('Enter the source files(ex:list.txt): ','s');
        lymph.dst=input('Enter the output file name(ex:/060112A/all): ','s');
        %Added on 03/20/17
        lymph.resol=input('Cut the spatial resolution half ?(y or n): ','s');
        
        lymph.src_loc=textread(lymph.src,'%s');
        lymph.dst_loc=lymph.dst;
        lymphatics_merge_image(lymph);
    case {7}
        lymph.src=input('Enter the source files(ex:list.txt): ','s');
        lymph.bas=input('Enter the baseline image file (ex:list.txt): ','s');
        lymph.src_loc=textread(lymph.src,'%s');
        lymph.bas_loc=textread(lymph.bas,'%s');
        lymphatics_percent_image(lymph);
    case {8}
        lymph.src=input('Enter the source file to calculate coreg (ex:src.img): ','s');
        lymph.modal=input('1.Adult_Rat 2.Baby_Rat 3.Prone 4. lateral 5. NIDA supine (Enter number): ');
        %lymph.ref=input('Enter the reference file  (ex:img_ref.img): ','s');
        
        if lymph.modal==1
        lymph.ref='/usr/local/matcodes/spm12_batch/Lymphatics_Prj/070912_template/pbase_snrv_Gadospin_062912A_E40.img';
        end
        if lymph.modal==2
        lymph.ref='/usr/local/matcodes/spm12_batch/Lymphatics_Prj/032613_baby_template/rrpbase_snrv_Magnevist_23_E65_032613.img';
        end
        if lymph.modal==3
        lymph.ref='/usr/local/matcodes/spm12_batch/Lymphatics_Prj/071714_Prone_template/rpbase_snrv_Magnevist_40_water_Cistern_prone_E50_071714A.img';
        end        
        if lymph.modal==4
        lymph.ref='/usr/local/matcodes/spm12_batch/Lymphatics_Prj/082014_Lateral_template/rpbase_snrv_Magnevist_40_water_Cistern_lateral_E71_082014A.img';
        end
        if lymph.modal==5
        lymph.ref='/usr/local/matcodes/spm12_batch/Lymphatics_Prj/072815_supine_template/rpbase_snrv_Magnevist_40_water_cistern_supine_NIDA_E55_072115B.nii';
        end
        
        lymph.oth=input('Enter the other file to apply the coregistration paramter (ex:other.txt): ','s');
        
        [lymph.dst lymph.tmp]=fileparts(lymph.src);
        
        lymph.oth_loc=textread(lymph.oth,'%s');
        lymphatics_realign_sagittal_image(lymph);
    case {9}  
        lymphatics_warp_batch;
    otherwise
        disp('Unknown option.')
end



