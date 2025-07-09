function lymphatics_realign_sagittal_image(lymph);

cd(lymph.dst);

spm_batch_src_dir='/usr/local/matcodes/spm12_batch/Lymphatics_Prj';
spm_batch_fname='lymphatics_spm_realign_sagittal_images.mat';

load([spm_batch_src_dir,'/',spm_batch_fname]);

matlabbatch{1}.spm.spatial.coreg.estwrite.ref{:}={};
matlabbatch{1}.spm.spatial.coreg.estwrite.source{:}={};
matlabbatch{1}.spm.spatial.coreg.estwrite.other{:}={};

matlabbatch{1}.spm.spatial.coreg.estwrite.ref{1}=[char(lymph.ref),',1'];
matlabbatch{1}.spm.spatial.coreg.estwrite.source{1}=[char(lymph.src),',1'];

for FN=1:size(lymph.oth_loc,1)
    matlabbatch{1}.spm.spatial.coreg.estwrite.other{FN,1}=[char(lymph.oth_loc(FN)),',1'];
end

%Run SPM
save([lymph.dst,'/',spm_batch_fname],'matlabbatch');
spm('defaults','pet');
spm_jobman('initcfg');
[job_id, idlist] = cfg_util('initjob',[lymph.dst,'/',spm_batch_fname]);
cfg_util('run',job_id);
cfg_util('deljob',job_id);
spm quit;

fprintf('Checkign affine.\n');

v_ref=spm_vol([char(lymph.ref)]);
img_ref=spm_read_vols(v_ref);

[fdir fname]=fileparts([char(lymph.src)]);
v_src=spm_vol([fdir,'/r',fname,'.nii']);
img_src=spm_read_vols(v_src);

for FN=1:size(lymph.oth_loc,1)
    clear CMD fname v_oth fdir img_oth img
    [fdir fname]=fileparts(char(lymph.oth_loc(FN)));
    v_oth=spm_vol([fdir,'/r',fname,'.nii']);
    
    if sum(sum((v_oth.mat+eps)./(v_ref.mat+eps))) ~= 16
        fprintf('Affine mismatched. Deleting r%s \n',fname);
        CMD=['rm ',fdir,'/r',fname,'{.nii}'];
        [status,result]=system(CMD);
    end
end

SL1=v_oth.dim(1)/2;
SL2=v_oth.dim(2)/2;
SL3=v_oth.dim(3)/2;

set(gcf,'Position',[1000 500 500 300])
subplot(2,3,1);imagesc(squeeze(img_ref(SL1,:,:)));colormap('gray');title('Reference');
subplot(2,3,2);imagesc(squeeze(img_ref(:,SL2,:)));colormap('gray');title('Reference');
subplot(2,3,3);imagesc(squeeze(img_ref(:,:,SL3)));colormap('gray');title('Reference');

subplot(2,3,4);imagesc(squeeze(img_src(SL1,:,:)));colormap('gray');title(['Source']);
subplot(2,3,5);imagesc(squeeze(img_src(:,SL2,:)));colormap('gray');title(['Source']);
subplot(2,3,6);imagesc(squeeze(img_src(:,:,SL3)));colormap('gray');title(['Source']);

%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3 6]);
print('-dpng',[lymph.dst,'/',fname,'.png']);

close all;

end


