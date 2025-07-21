function lymphatics_percent_image(lymph)



%Create mean image based on the baselines 
for FN=1:size(lymph.bas_loc,1)
    clear v img;

    v=spm_vol(char(lymph.bas_loc(FN)));
    img=spm_read_vols(v);

    if FN==1
        img_bas=img;
    else
        img_bas=img_bas+img;
    end

end


img_bas=img_bas/size(lymph.bas_loc,1);

img_sum=zeros(size(img_bas));


for FN=1:size(lymph.src_loc,1)
    clear v img fdir fname;
    v=spm_vol(char(lymph.src_loc(FN)));
    img=spm_read_vols(v);

    %Normalize it.
    img=(img-img_bas+eps)./(img_bas+eps)*100;
    [fdir fname]=fileparts(char(lymph.src_loc(FN)));
    v.fname=[fdir,'/p',fname,'.nii'];
    spm_write_vol(v,img);
    fprintf('%s done \n',fname);

    %Sum all images
    img_sum=img_sum+img;
end


v.fname=[fdir,'/pbase_',fname,'.nii'];
spm_write_vol(v,img_bas);
v.fname=[fdir,'/psum_',fname,'.nii'];
spm_write_vol(v,img_sum);


end