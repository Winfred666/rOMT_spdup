function lymphatics_realign_image(lymph)

% old_dir = pwd;

spm_batch_src_dir='/home/xym/Desktop/MRI_ROMT/rOMT_spdup-main/data/ours/config';
spm_batch_fname='realign_batch_config.mat';

S = load([spm_batch_src_dir,'/',spm_batch_fname], 'matlabbatch');
matlabbatch = S.matlabbatch;
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


fprintf("checkout outcome in %s\n", lymph.dst);
% cd(old_dir);

end



