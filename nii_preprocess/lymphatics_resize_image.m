function lymphatics_resize_image(lymph)

for FN=1:size(lymph.src_loc,1)
    clear v img fname v_mat;
    fname=char(lymph.src_loc(FN));
    v=spm_vol(fname);
    img=spm_read_vols(v);
    v_mat=spm_imatrix(v.mat);
    
    %added on 4/8/18
    if lymph.species==1  %rat
        v_mat(7:9)=v_mat(7:9)*10; %Multiply the resolution by 10
    elseif lymph.species==2 %mouse
        v_mat(7:9)=v_mat(7:9)*20; %Multiply the resolution by 20
    end
    
    
    v.mat=spm_matrix(v_mat);

    v.fname=char(lymph.dst_loc(FN));%Assign destination file name
    spm_write_vol(v,img);
    [tmp_dir tmp_fname]=fileparts(fname);

    fprintf('%s done \n',tmp_fname);
end


end