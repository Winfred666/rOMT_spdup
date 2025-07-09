function lymphatics_warp_apply(lymph)
clear v spat_dim;

cd(lymph.dst);

spm_batch_src_dir='/usr/local/matcodes/hedok/spm8_batch/Lymphatics_Prj';
spm_batch_fname='lymphatics_spm_warp_apply.mat';

load([spm_batch_src_dir,'/',spm_batch_fname]);


matlabbatch{1}.spm.spatial.normalise.write.subj.matname={lymph.spmnorm};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample{1}=[lymph.pbas,',1'];
matlabbatch{1}.spm.spatial.normalise.write.subj.resample{2}=[lymph.psum,',1'];

v=spm_vol(char(lymph.mbpbas));
spat_dim=spm_imatrix(v.mat);
matlabbatch{1}.spm.spatial.normalise.write.roptions.vox=abs(spat_dim(7:9));


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