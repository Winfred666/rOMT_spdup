function lymphatics_sum_image(lymph)

img_all=[];

for FN=1:size(lymph.src_loc,1)
    clear v img fname v_mat;
    fname=char(lymph.src_loc(FN));
    v=spm_vol(fname);
    img=spm_read_vols(v);

    if FN==1
        img_all=img;
    else
        img_all=img_all+img;
    end
    [tmp_dir tmp_fname]=fileparts(fname);
    fprintf('Adding %s \n',tmp_fname);
end

img_all=img_all/size(lymph.src_loc,1); %Take the mean

img_all(isnan(img_all))=0;

v.fname=char(lymph.dst_loc);%Assign destination file name
spm_write_vol(v,img_all);

[tmp.dir tmp.fname]=fileparts(lymph.dst_loc);

save([tmp.dir,'/lymphatics_sum_img'],'lymph');


end