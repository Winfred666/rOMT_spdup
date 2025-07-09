function lymphatics_warp_normalize(lymph)
clear v spat_dim;

cd(lymph.dst);

spm_batch_src_dir='/usr/local/matcodes/hedok/spm8_batch/Lymphatics_Prj';
spm_batch_fname='lymphatics_spm_warp_normalize.mat';
lymph.template='/usr/local/matcodes/hedok/spm8_batch/Lymphatics_Prj/normalize_template/mbpbase_snrv_Magnevist_23_E65_061313.nii';

load([spm_batch_src_dir,'/',spm_batch_fname]);
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.source={};
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample={};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.template={};

v=spm_vol(char(lymph.src));
spat_dim=spm_imatrix(v.mat);

matlabbatch{1}.spm.spatial.normalise.estwrite.subj.source={[char(lymph.src),',1']};
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample={[char(lymph.src),',1']};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.template={[char(lymph.template),',1']};
matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.vox=abs(spat_dim(7:9));


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