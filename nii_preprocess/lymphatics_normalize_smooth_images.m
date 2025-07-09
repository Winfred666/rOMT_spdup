function lymphatics_normalize_images(lymph)

%Delete the destination directory if it exists.
if exist(lymph.dst,'dir')
    CMD=['rm -rf ',lymph.dst];
end

CMD=['mkdir ',lymph.dst];
[status,result]=system(CMD);

cd(lymph.dst);

clear v_msk img_msk
v_msk=spm_vol(lymph.msk);
img_msk=spm_read_vols(v_msk);
img_msk=img_msk>0;

for FN=1:size(lymph.src_loc,1)
    clear fname v img img_scale mean_intensity img_norm tmp fname_out sfname_out SG;
    
    fname=char(lymph.src_loc(FN));
    v=spm_vol(fname);
    img=spm_read_vols(v);
    
    if sum(sum((v_msk.mat(1:3,1:3)+eps)./(v.mat(1:3,1:3)+eps)))==9 %Check and make sure affine in mask is the same as all images
        
        img_scale=img.*img_msk;
        mean_intensity=mean(img_scale(img_scale>0));
        img_norm=img/mean_intensity*1000; %normalize
        
        %Normalize image
        [tmp.dir tmp.fname]=fileparts(fname);
        fname_out=[lymph.dst,'/n',tmp.fname,'.nii'];
        v.fname=fname_out;
        spm_write_vol(v,img_norm);
        
        %Smooth the image
        sfname_out=[lymph.dst,'/sn',tmp.fname,'.nii'];
        SG=[lymph.smooth lymph.smooth lymph.smooth];
        spm_smooth(fname_out,sfname_out,SG)
        
        %Added on 03/20/17
        if strcmp(lymph.print,'y')
            clf;
            set(gcf,'Position',[1000 500 500 1000])
            clear img2d
            img2d=squeeze(img(:,lymph.SL,:));
            img2d_norm=squeeze(img_norm(:,lymph.SL,:));
            
            maxI1=max(img2d(img2d>0));minI1=maxI1*0.1;
            minI2=150;maxI2=2500;
            
            tmp.fname=strrep(tmp.fname,'_','-');
            subplot(2,1,1);imagesc(img2d);caxis([minI1 maxI1]);title(tmp.fname);
            subplot(2,1,2);imagesc(img2d_norm);caxis([minI2 maxI2]);title(['Normalized ',tmp.fname]);
            
            pause(1)
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3 6]);
            print('-dpng',[lymph.dst,'/',tmp.fname,'.png']);
        end
        fprintf('%s done. \n',fname);
        
        
    else
        pause
        fprintf('%s affine transformations do not match with the mask in normalization \n',fname);
        
    end
    
end

    save([lymph.dst,'/Normalization_info'],'lymph');


end