function lymphatics_warp_relice(lymph);

cd(lymph.dst);

spm_batch_src_dir='/usr/local/matcodes/hedok/spm8_batch/Lymphatics_Prj';
spm_batch_fname='lymphatics_spm_warp_reslice.mat';

load([spm_batch_src_dir,'/',spm_batch_fname]);
matlabbatch{1}.spm.spatial.coreg.write.ref={};
matlabbatch{1}.spm.spatial.coreg.write.source={};

matlabbatch{1}.spm.spatial.coreg.write.ref=lymph.ref;
matlabbatch{1}.spm.spatial.coreg.write.source=lymph.src;

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