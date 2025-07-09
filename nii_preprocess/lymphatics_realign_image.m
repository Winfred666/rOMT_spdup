function lymphatics_realign_image(lymph);

cd(lymph.dst);

spm_batch_src_dir='/usr/local/hedok/matcodes/spm12_batch/Lymphatics_Prj';
spm_batch_fname='lymphatics_spm_realign_images.mat';

load([spm_batch_src_dir,'/',spm_batch_fname]);
matlabbatch{1}.spm.spatial.realign.estwrite.data{:}={};

for FN=1:size(lymph.src_loc,1)
    matlabbatch{1}.spm.spatial.realign.estwrite.data{1}{FN,1}=[char(lymph.src_loc(FN)),',1'];
end

%Run SPM
save([lymph.dst,'/',spm_batch_fname],'matlabbatch');
spm('defaults','pet');
spm_jobman('initcfg');
[job_id, idlist] = cfg_util('initjob',[lymph.dst,'/',spm_batch_fname]);
cfg_util('run',job_id);
cfg_util('deljob',job_id);
spm quit;
%

end



