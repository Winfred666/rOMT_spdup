clc;
clear all;


fprintf('1. Extract brain using the coregisterd mask \n' );
fprintf('2. Affine transform B1 corrected images \n' );
fprintf('3. Apply transform pbase and psum images \n' );
fprintf('4. Take mean pbase psum images \n' );

run_N=input('Enter the numer:');

switch lower(run_N)
    case {1}
        src_loc=input('Enter the list of mask file names: ','s');
        ref_loc=input('Enter the list of pbase file names: ','s');

        lymph.src_loc=textread(src_loc,'%s');
        lymph.ref_loc=textread(ref_loc,'%s');

        for FN=1:size(lymph.src_loc,1)
            clear msk_dir msk_fname v_msk img_msk v img file_dir file_name;

            %Reslice the mask
            [lymph.dst tmp]=fileparts(char(lymph.ref_loc(FN)));
            lymph.ref=lymph.ref_loc(FN);
            lymph.src=lymph.src_loc(FN);
            lymphatics_warp_reslice(lymph);

            %Binarize the mask
            [msk_dir msk_fname]=fileparts(char(lymph.src));
            rmsk_fname=[msk_dir,'/r',msk_fname,'.img'];
            v_msk=spm_vol(rmsk_fname);
            img_msk=spm_read_vols(v_msk);
            img_msk(img_msk>0)=1;

            v=spm_vol(char(lymph.ref));
            img=spm_read_vols(v);
            img=img.*img_msk;
            [file_dir file_name]=fileparts(v.fname);
            v.fname=[file_dir,'/b',file_name,'.img'];
            spm_write_vol(v,img);

            FN
        end
    case {2}
        src_loc=input('Enter the list of B1 corr and brain file names: ','s');
 
        lymph.src_loc=textread(src_loc,'%s');
        
        for FN=1:size(lymph.src_loc,1)
            lymph.src=lymph.src_loc(FN);
            [lymph.dst tmp]=fileparts(char(lymph.src_loc(FN)));
            lymphatics_warp_normalize(lymph);
        end
        FN
    case {3}
        pbas_loc=input('Enter the list of pbase file names: ','s');
        psum_loc=input('Enter the list of psum file names: ','s');
        mbpbas_loc=input('Enter the list of mbpbase file names: ','s');
        spmnorm_loc=input('Enter the list of normalization parameters: ','s');
        
        
        lymph.pbas_loc=textread(pbas_loc,'%s');
        lymph.psum_loc=textread(psum_loc,'%s');
        lymph.mbpbas_loc=textread(mbpbas_loc,'%s');
        lymph.spmnorm_loc=textread(spmnorm_loc,'%s');

        for FN=1:size(lymph.pbas_loc,1)
            clear v_pbas v_psum v_mbpbas pbas_dir pbas_fname psum_dir psum_fname;
            
            v_pbas=spm_vol(char(lymph.pbas_loc(FN)));
            v_psum=spm_vol(char(lymph.psum_loc(FN)));
            v_mbpbas=spm_vol(char(lymph.mbpbas_loc(FN)));
            
            lymph.mbpbas=char(lymph.mbpbas_loc(FN));
            lymph.spmnorm=char(lymph.spmnorm_loc(FN));
            
            img_pbas=spm_read_vols(v_pbas);
            img_psum=spm_read_vols(v_psum);
            img_mbpbas=spm_read_vols(v_mbpbas);
            
            [pbas_dir,pbas_fname]=fileparts(v_pbas.fname);
            [psum_dir,psum_fname]=fileparts(v_psum.fname);
             
            %Create new pbas
            img_mbpbas(:)=0;
            img_mbpbas=img_pbas;
            v_mbpbas.fname=[pbas_dir,'/Aff_',pbas_fname,'.nii'];
            spm_write_vol(v_mbpbas,img_mbpbas);
            lymph.pbas=v_mbpbas.fname;
            
            %Create new psum
            img_mbpbas(:)=0;
            img_mbpbas=img_psum;
            v_mbpbas.fname=[psum_dir,'/Aff_',psum_fname,'.nii'];
            spm_write_vol(v_mbpbas,img_mbpbas);
            lymph.psum=v_mbpbas.fname;
            
            [lymph.dst tmp]=fileparts(char(lymph.mbpbas_loc(FN)));
            lymph.mbpbas=char(lymph.mbpbas_loc(FN));
            
            lymphatics_warp_apply(lymph)
            FN
        end
    case {4}
        pbas_loc=input('Enter the list of warped pbase file names: ','s');
        ave_pbas_fname=input('Enter the output file_name for pbase: ','s');
        ave_psum_fname=input('Enter the output file_name for psum: ','s');
        
        lymph.pbas_loc=textread(pbas_loc,'%s');
        
        file_N=size(lymph.pbas_loc,1);
        
        err=0;
        img_tot_pbas=zeros([256,256,256]);
        img_tot_psum=zeros([256,256,256]);
        
        
        for FN=1:size(lymph.pbas_loc,1)
            clear v_pbas v_psum img_pbas img_psum;
            
            v_pbas=spm_vol(char(lymph.pbas_loc(FN)));
            img_pbas=spm_read_vols(v_pbas);

            fname_psum=strrep(char(lymph.pbas_loc(FN)),'base','sum');

            v_psum=spm_vol(fname_psum);
            img_psum=spm_read_vols(v_psum);
                        
            if FN==1
                dir_cos=v_pbas.mat;
            end

            x_dim=size(img_pbas,1);
            y_dim=size(img_pbas,2);
            z_dim=size(img_pbas,3);

            %Some files are 256x255x255 and other are 256x256x256 so this is necessary
            img_tot_pbas([1:x_dim],[1:y_dim],[1:z_dim])=img_tot_pbas([1:x_dim],[1:y_dim],[1:z_dim])+img_pbas;
            img_tot_psum([1:x_dim],[1:y_dim],[1:z_dim])=img_tot_psum([1:x_dim],[1:y_dim],[1:z_dim])+img_psum ;
           
            
            if round(sum(sum((dir_cos+eps)./(v_pbas.mat+eps)))/16)~=1
                fprintf('Direciton of cosine does not match \n');
                err=1;
            end
            
        end

        if err==0
            img_tot_pbas=img_tot_pbas/file_N;
            img_tot_psum=img_tot_psum/file_N;
            
            v_pbas.dim=size(img_tot_pbas);
            v_psum.dim=size(img_tot_psum);

            v_pbas.fname=ave_pbas_fname;
            spm_write_vol(v_pbas,img_tot_pbas);

            v_psum.fname=ave_psum_fname;
            spm_write_vol(v_psum,img_tot_psum);
        end

        
    otherwise


        disp('Unknown option.')
end


