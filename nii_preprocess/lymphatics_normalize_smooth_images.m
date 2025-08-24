function lymphatics_normalize_smooth_images(lymph)

%Delete the destination directory if it exists.
if exist(lymph.dst,'dir')
    CMD=['rm -rf ',lymph.dst];
    system(CMD);
end

CMD=['mkdir ',lymph.dst];
[status,result]=system(CMD);

% cd(lymph.dst);

clear v_msk img_msk
v_msk=spm_vol(lymph.msk);
img_msk=spm_read_vols(v_msk);
img_msk=img_msk>0.01;

for FN=1:size(lymph.src_loc,1)
    clear fname v img img_scale mean_intensity img_norm tmp fname_out sfname_out SG;
    % get the folder of src image
    fname=char(lymph.src_loc(FN));
    v=spm_vol(fname);
    img=spm_read_vols(v);
    
    % resample and align header if fail to align!!
    affine_ratio = (v_msk.mat(1:3,1:3)+eps)./(v.mat(1:3,1:3)+eps);
    if all(abs(abs(affine_ratio(:)) - 1) < 1e-3)
        img_scale=img.*img_msk;
        mean_intensity=mean(img_scale(img_scale>0));
        img_norm=img/mean_intensity*1000; %normalize
        
        %Normalize image and save
        [tmp.dir, tmp.fname]=fileparts(fname);
        
        fname_out=[tmp.dir,'/n',tmp.fname,'.nii'];
        v.fname=fname_out;
        if exist(lymph.dst,'file')  % if the file already exists, delete it
            system(['rm -f ', char(fname_out)]);
        end

        spm_write_vol(v,img_norm);
        
        %Smooth the image and save again
        sfname_out=[tmp.dir,'/sn',tmp.fname,'.nii'];
        SG=[lymph.smooth lymph.smooth lymph.smooth];
        spm_smooth(fname_out,sfname_out,SG) % still save into src.

        system(['rm -f ', char(fname_out)]); % only leave sfname_out

        %Added on 03/20/17
        if strcmp(lymph.print,'y')
            fig = figure('Visible', 'off'); % Create an invisible figure
            set(fig,'Position',[1000 500 500 1000])
            clear img2d
            img2d=squeeze(img(:,lymph.SL,:));
            img2d_norm=squeeze(img_norm(:,lymph.SL,:));
            
            maxI1=max(img2d(img2d>0));minI1=maxI1*0.1;
            minI2=150;maxI2=2500;
            
            tmp.fname=strrep(tmp.fname,'_','-');
            subplot(2,1,1);imagesc(img2d);caxis([minI1 maxI1]);title(tmp.fname);
            subplot(2,1,2);imagesc(img2d_norm);caxis([minI2 maxI2]);title(['Normalized ',tmp.fname]);
            
            set(fig,'PaperUnits','inches','PaperPosition',[0 0 3 6]);
            print(fig, '-dpng', [lymph.dst,'/',tmp.fname,'.png']);
            close(fig); % Close the figure to free up memory
        end
        fprintf('%s done. \n',fname);
        
        
    else
        % pause
        fprintf('%s affine transformations do not match with the mask in normalization \n',fname);
        
    end
    
end

    save([lymph.dst,'/Normalization_info'],'lymph');


end