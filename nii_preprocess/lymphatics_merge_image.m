function lymphatics_merge_image(lymph)

[fout_dir fout_name]=fileparts(lymph.dst);


%Reduce the resolution by half so it's faster       

for FN=1:size(lymph.src_loc,1)
    clear CMD fname fdir img2 img
    [fdir fname]=fileparts(char(lymph.src_loc(FN)));
    v=spm_vol([fdir,'/',fname,'.nii']);
    img=spm_read_vols(v);
    
    if strcmp(lymph.resol,'y')
        if FN==1
            fprintf('Reducing resolution by half.\n');
        end        
        img2=img(1:2:end,1:2:end,1:2:end);
        v.dim=v.dim/2;
        v.fname=[fdir,'/tmp_',fname,'.nii'];
        spm_write_vol(v,img2);        
    else
        img2=img;
        v.dim=v.dim;
        v.fname=[fdir,'/tmp_',fname,'.nii'];
        spm_write_vol(v,img2);                
    end
end

for FN=1:size(lymph.src_loc,1)
    clear CMD fname fdir
    [fdir fname]=fileparts(char(lymph.src_loc(FN)));
    
    if FN==1
        CMD=['cp ',fdir,'/tmp_',fname,'.nii ',fout_dir,'/',fout_name,'.nii'];
        [status,result]=system(CMD);
        CMD=['fslchfiletype NIFTI_GZ',' ',fout_dir,'/',fout_name];
        [status,result]=system(CMD);
    else
        CMD=['fslmerge -t ',fout_dir,'/',fout_name,' ',fout_dir,'/',fout_name,' ',fdir,'/tmp_',fname];
        [status,result]=system(CMD);

    end
    fprintf('%s done!\n',fname);
    
    CMD=['rm ',fdir,'/tmp_',fname,'.nii'];
    [status,result]=system(CMD);
    
end



end