
# 2025 调整版本运行方式

## 1. 数据预处理

脚本位于 `lymphatics_batch.m`，使用 spm_batch 做处理，必要步骤为 
3：remove head motion， 
5：normalize and smooth each frame
7: get percentage on baseline；

### 1.1 准备 DCE-MRI 数据文件夹

保存好 baseline，src，里面只有 `.nii.gz` 后缀文件合法
```matlab
lymph.bas_loc_dir = '/data/xym/DEX_MRI/ISO/ISO_52/DCE_nii_baseline';
lymph.src_loc_dir = '/data/xym/DEX_MRI/ISO/ISO_52/DCE_nii_data';
```

### 1.2 准备 Template.nii

空间对齐，并且 nii 头文件变换矩阵也对齐（否则报错）的全脑掩膜，此处预处理是用于 确定 normalize 的平均值取值范围；之后 ROMT 也要用该掩膜。

用 C57Bl6 ，通过 3D Slicer General Register + Resample ，在文件中设置：

```matlab
lymph.msk = '/data/xym/DEX_MRI/ISO/Template_C57Bl6_n30_brain_ISO_52.nii';
```

### 1.3 运行 lymphatics_batch.m

所有参数调整都在 `lymphatics_batch.m` 文件中；可调整的参数有：

```matlab
lymph.force_mass_conservation   % 是否强制质量守恒，每一帧
lymph.smooth = 0.125 * 2;       % 三维平滑滤波 掩膜/核 的大小，为体素单位 * 核的体素宽度
```

之后运行：

```sh
matlab -nodisplay -nosplash -r "run('lymphatics_batch.m')"
```

程序运行后，会在与源文件相同的 `lymph.src` 目录下得到 `psnr_` 开头的预处理文件。


### 1.4 准备 DTI 数据

较好的处理流程：


在 dicom 文件夹中：

1. 

```bash
mkdir -p out_dwi
```

2. 
```bash
dcm2niix -o out_dwi -f dwi_%p_%s -b y -z y ./E*_DTI_sat_30
```
利用原始扫描数据，输出扫描空间下的 DWI 系数场， b-value 和 b-vector

3. 

```bash
fslval out_dwi/dwi_DTI_EPI_30dir_sat_*.nii.gz dim4      # number of volumes
wc -w out_dwi/dwi_DTI_EPI_30dir_sat_*.bval              # must match dim4
```
快速检查输出文件所含 b 值个数

4. 

```bash
cd out_dwi;dwigradcheck dwi_DTI_EPI_30dir_sat_*.nii.gz -fslgrad dwi_DTI_EPI_30dir_sat_*.bvec dwi_DTI_EPI_30dir_sat_*.bval
```

检查 bval 和 bvec gradients，若不成功需要用 dwigradcorrect

5. 
```bash
fslroi dwi_*.nii.gz b0 0 1
```
使用 b0 做大脑mask；用 cat dwi_DTI_EPI_30dir_sat_*.bval 查看 diffuse 最小权重作为 b0 ，通常是第一个，权重保存为 “Template_T2_DTI_brain.nii.gz”；
bet b0 b0_brain -m  # 去颅骨做 mask，难以得到准确 mask，从C57Bl6_T2手动配准获取Template_T2_DTI_brain.nii.gz

6. 

```bash
T1ref="rT1_FLASH_3D_0020.nii"
b0_image="b0.nii.gz"

# flirt -in Template_T2_DTI_brain.nii.gz \
#      -ref "${b0_image}" \
#      -applyxfm -usesqform \
#      -out mask_in_dwi.nii.gz
```

3D slicer 的 Linear / Bspline resample 精度不高，但可以使用 WindowedSinc 确保对齐 header ，精度较高，不用 flirt


7. 
```bash
dtifit -k dwi_DTI_EPI_30dir_sat_*.nii.gz \
       -o dti \
       -m mask_in_dwi.nii.gz \
       -r dwi_DTI_EPI_30dir_sat_*.bvec \
       -b dwi_DTI_EPI_30dir_sat_*.bval \
       --save_tensor
```
利用每一个方向的 bval 和 bvec 计算 DTI ，得到 DTI 的特征方向 V1,V2,V3，特征值 L1,L2,L3,以及6D对称矩阵场 dti_tensor.nii.gz；

8. 
```bash
ImageMath 3 dti_ants.nii.gz FSLTensorToITK dti_tensor.nii.gz
```
一定要注意FSL 得到的每个 voxel 6D tensor 在 FSL 中是上三角排列：Dxx Dxy Dxz Dyy Dyz Dzz，必须转为下三角，适合 ANTs/ITK Format: Dxx Dxy Dyy Dxz Dyz Dzz；同时 header 是 NIFTI-1 format；数据格式变为 5D x-y-z-1-6
参考：https://github.com/ANTsX/ANTs/wiki/Importing-diffusion-tensor-data-from-other-software

9. 
```bash
mkdir -p dti_aligned
antsRegistration --dimensionality 3 \
    --output [dwi_to_t1_, dti_aligned/b0.nii.gz] \
    --winsorize-image-intensities [0.005,0.95] \
    --initial-moving-transform "[${T1ref},${b0_image},1]" \
    --transform Affine[0.1] \
    --metric MI[${T1ref},${b0_image},1,32] \
    --convergence [100x70x50,1e-6,10] \
    --shrink-factors 4x2x1 \
    --smoothing-sigmas 2x1x0vox
```
类似 flirt 但参数更细，直接得到 dwi_to_t1_0GenericAffine.mat


10. 
```bash
antsApplyTransforms -d 3 \
                    -e 2 \
                    -i dti_ants.nii.gz \
                    -r "$T1ref" \
                    -t dwi_to_t1_0GenericAffine.mat \
                    -o dti_aligned/dti_tensor.nii.gz \
                    -f 0.0007 \
                    --verbose
```
注意对于 tensor ，默认的就是 Log-Euclidean interpolation，不需要再调整 -n or --interpolation
而 -f 0.0007 设置背景的扩散系数默认值，防止 log(0)，非常有用

11. 这时会生成一个 `dti_aligned/dti_tensor.nii.gz` 文件。可回到 matlab 中，设置 `base_path`，运行`dti_preprocess/genDmatFromDtifit.m`：

```matlab
addpath('./dti_preprocess');
addpath('./utilities');
base_path = '/data/xym/DTI_data/KX 078/dicom/out_dwi/dti_aligned';
genDmatFromDtifit
```

之后便得到 `dti_aligned/dti_tensor.mat` ，可作为算法配置，同时也可用代码检查 DTI 的 MD 和 FA 值是否正确。

```matlab
dti_sanity_check
```

## 2. 运行 ROMT 优化

### 2.1 参数调整

ROMT 优化过程的参数调整，可在 `configs` 文件夹中，复制一份已有的 json 文件；在 `set_config_CAA.m` 会将里面的所有参数载入。

设置 `only_post_processing = 0` 以执行 ROMT 优化。以下是优化过程的具体参数：

关于 DTI 的参数：

```matlab
dti_path        % 填入之前 dti_aligned/dti_tensor.mat 的绝对路径
dti_enhanced    % 对 DTI 信号的增强倍数，若可视化DTI 时，网格中大部分 MD 值小于 2e-3 时，可用此强化扩散效应。
```

基础参数：

```matlab
data_template   % 之前预处理得到的 psnr_ 序列
ROI_msk_path    % 全脑分析中，使用之前预处理使用的 Template.nii 即可

do_ROI_msk，ROI_msk_threshold  % 是否启用上述 mask，值大于 threshold 的体素考虑在 mask 内。 


x_range, y_range, z_range % 与数据保持一致，如nii数据为x轴方向分辨率为 128，那么设置 x_range 为 [1,128]；

do_resize       % 设置为 1 时会使用 `size_factor` 下采样数据，显著提高运行速度

smooth          % 在空间维度上再次平滑密度场，设置为 dt 的几倍表示平滑范围为多少格；

dilate          % ROI_msk 膨胀 dilate 个体素后才是优化的最终控制域，防止边界不平滑处流动受限，

first_time      % 从该数据帧开始 优化/插值
time_jump       % 设置为 2，表示每次 优化/插值 的初始密度场为 sprintf("%04d.nii",ti)，目标密度场为 sprintf("%04d.nii",ti+2) 
last_time       % 该数据帧是最后一次 优化/插值 的初始密度场（因此需要 last_time + time_jump 帧）

sigma           % 关键：经验扩散系数
dx              % 每个体素的长度，单位毫米，只是用于将 DTI 从真实物理单位换算成格子单位。
dt              % 关键：插值帧之间（不是原始数据帧之间）的最小时间步长
nt              % 关键：插值帧数，一般设置 nt > 1 以增加自由度，一次优化过程 = 一次产生 nt 个中间帧的插值过程

% 如果MRI间隔为4分钟，确保nt*dt=4以获得分钟的时间分辨率

gamma,beta  % loss 权重，不用动
reinitR,reInitializeU  % 重置而非使用上一次优化得到的密度/速度场。并行则都设置为 1， 30 min 可跑完；串行都设置 0， 5 hour 才跑完但得到速度更快更平滑；

niter_pcg % 每次更新线性方程求解/优化更新的迭代次数，一个递增数组，当低次迭代失败时，会尝试更高迭代次数，一般只填入一个数字即可，如 [500]。
% 对于 20^3 网格（下采样0.2）设置为 [80]；对于 60^3 网格（下采样 0.5） 设置为 [500]；100^3 网格设置为 [700] 较合适。虽然会非常慢，需要跑一到两天，但求解精准。relres 相对残差下降到 0.1 以下才合格。 

maxUiter  % 建议将 niter_pcg 调大，这样 maxUiter 设置 20 即可收敛
```

### 2.2 运行

若运行 ISO 配置的优化过程，则执行：

```sh
matlab -nodisplay -nosplash -r "config_path='ours_ISO';run('driver_CAA.m')"
```

其中 config_path 与配置的 json 文件名相同。

## 3. 运行 GLAD 后处理

### 3.1 参数调整

若之前执行过同配置的 ROMT 优化，则设置 `only_post_processing = 1` ，以从保存的 `u0_ours_ISO...mat` 中恢复速度场，只进行后处理。

```matlab
exclude_frames  % 在后处理中过滤异常速度场，比如 4 表示设置 0004-0005 的速度为 0，不移动轨迹线或用于计算平均速度

density_percent_thres   = 16; % 过滤全过程所有时刻中，最大相对密度低于此阈值的体素的速度场，原同 exclude_frames

sp_thresh                     % 过滤全过程所有时刻中，最大相对密度低于此阈值的体素，不在这些体素生成 pathline 起始点

sl_tol                        % 只保留长度多于 sl_tol 格的 pathline

GLAD_visualize_mask: {
    enable:True,
    path: "...",
    threshold: "0.01",
}         % GLAD 可视化 pathline, speed map, flux vector 等等图像时使用的背景板。只会保留从这个 mask 出发的 pathline，可用于看清某个具体位置的 pathline 情况。


compartment_mask_list: [
    {
        name: "White matter",
        path: "...",
        threshold:"0.01",
    }
]       % 分别对多个 mask 中的脑区计算该脑区中的平均速度、平均 Pe 值、v-flux（覆盖体积） 值

GLAD_timestep_factor        % pathline 位移延长，调为 1 pathline 才符合真实。
```


其余 `speedmap_slice` `strid` 等，是最后可视化参数调整，如果不对，完全可以跑完 GLAD 后，将 `driver_CAA.m` 内部的绘图代码，复制到命令行中快速重跑可视化。

此外在 `paramInitGLADpar.m` 中，有利用 `smoothn.m` 对速度场和 pathline 进行平滑的参数 `glaSvt`， `glaSvs`， `glasmpTol`，速度和 pathline 不平滑可适当调大些。

### 3.2 运行

设置 `only_post_processing = 1` 后，重新跑：

```sh
matlab -nodisplay -nosplash -r "config_path='ours_ISO';run('driver_CAA.m')"
```




> 原版 README 见：
>
> [Introduction of rOMT_spdup](https://github.com/xinan-nancy-chen/rOMT_spdup)


## 4. 已有速度场做可视化

当 `only_post_processing` 为 1 时，默认会在 `output_dir` 中读取速度场并运行 GLAD ，也可以指定速度场路径。

要求速度场的格式是 `.mat` ，其中含有 N x 1 的向量 `u`。假设一共使用了 ft = (last_time - first_time)/time_jump 帧的密度场，则 N = (x_range * y_range * z_range) * 3 * nt * ft。

`u` 形如

```
{[u_111,... ,u_xyz, v_111,... ,v_xyz, w_111,... ,w_xyz] for t_1_1, ... , [u_111,... ,u_xyz, v_111,... ,v_xyz, w_111,... ,w_xyz] for t_1_nt}, ... , { ... [u_111,... ,u_xyz, v_111,... ,v_xyz, w_111,... ,w_xyz] for t_ft_nt}
```

注意密度场的配置如 `size_factor` ，要与其他方法的运行条件相同，将 `only_post_processing` 设置为 1，然后将 `steady_velocity_file_path` 设置为该速度场文件路径即可。


## 5. 注意事项

长度单位统一为 cell 即设 DCEMRI 体素为单位 1 ，不是真实单位；时间单位统一为 minute ，如果 MRI 记录间隔 4min，那么需要保证 nt * dt = 4 ；而载入的 DTI 默认单位为 mm^2/s，在计算之前也已经换算成 cell^2/min。需要自行在保存的 `test_results` 中将数据换算回真实单位；
