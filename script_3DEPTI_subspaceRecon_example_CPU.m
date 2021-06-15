% Example script of the low-rank subspace reconstruction for 3D-EPTI using BART on CPU
% Related to the following works
% 1) "3D Echo Planar Time-resolved Imaging (3D-EPTI) for ultrafast multi-parametric quantitative MRI",Fuyixue Wang and Zijing Dong et al., bioRxiv 2021.
% 2) "Echo planar time-resolved imaging (EPTI)", Fuyixue Wang et al.,MRM,2019
% 3) "EPTI with Subspace Reconstruction and Optimized Spatiotemporal Encoding", Zijing Dong et al.2020
% Fuyixue Wang <fwang18@mgh.harvard.edu>, Zijing Dong <zijingd@mit.edu> Feb/2021

clear;
close all;
addpath('Funcs');
bart_path = ''; % please add the BART path here, BART needs to be installed correctly
setenv('TOOLBOX_PATH', bart_path)
addpath(strcat(getenv('TOOLBOX_PATH'), '/matlab'));
%% Data Preparation for bart
dir_BART = ['tmp/'];
mkdir(dir_BART);
K = 12;
load('Data_forRecon.mat')
load('Generated_bases.mat')

Phi = U(:,1:K);
[npe,nz,~,~] = size(kdata);
ISize=FillOnesTo16([npe,nz,1,1,1,K]);
Phi=permute(Phi,[3 4 5 1 6 2]); 

[scaling] = Esti_Scaling(kdata);
fprintf('\nScaling: %f\n\n', scaling);

sens_FN=[dir_BART,'sens'];
writecfl(sens_FN,sens_map_tmp); % save sensitivity maps
kdata_FN=[dir_BART,'kdata'];
writecfl(kdata_FN,single(kdata)); % undersampled k-t data
Phi_use = repmat(Phi,[npe nz]);
Phi_FN=[dir_BART,'PhiUse'];       % subspace bases
writecfl(Phi_FN,Phi_use);
Phase_FN=[dir_BART,'Phase_T'];    % phase evolution maps obtained from calibration scan
writecfl(Phase_FN,Phase_T);
sample_mask_FN=[dir_BART,'sample_mask'];     % undersampling mask
writecfl(sample_mask_FN,single(mask_sample));
    
ScriptFN=[pwd filesep 'EPTI_CompsToSig.txt'];
WriteLinopToFile(ScriptFN,{'fmac 0 32','fmac 1 0','fmac 2 0','fftc 3','fmac 3 0'}); % (1) x Phi and sum along dim6 (2) x Phase (3) x Sens (4) FFT in dim1&2 (5) x Mask

%% Reconstruction
disp('Reconstruction Start');
tic;    
Rec=bart(['picsS -S -d 5 -i 100 -u 0.01 -C 4 -m -w ',num2str(scaling),' -b 10 -R L:3:3:0.0003 ' ScriptFN],ISize,kdata_FN,Phi_FN,Phase_FN,sens_FN,sample_mask_FN); % CPU LLR
toc;
res_a = Rec;
im_recon=squeeze(sum(Rec.*Phi,6));
%% Show results
nt_GE = 53;
nacq_GE = 20;
nt_SE = 43;
nacq_SE = 10;

im_recon = permute(im_recon,[2,1,3]);
im_recon_GE = im_recon(:,:,1:nt_GE*nacq_GE);
im_recon_SE = im_recon(:,:,nt_GE*nacq_GE+1:end);

im_recon_GE=reshape(im_recon_GE,[size(im_recon_GE,1),size(im_recon_GE,2),nt_GE,nacq_GE]);
im_recon_SE=reshape(im_recon_SE,[size(im_recon_SE,1),size(im_recon_SE,2),nt_SE,nacq_SE]);

T_to_show_GE=3:3:nacq_GE;
Echo_forPhase_GE = ceil(nt_GE/2)-10:ceil(nt_GE/2)+10;
T_to_show_SE=3:2:8;
Echo_forPhase_SE = ceil(nt_SE/2)-5:ceil(nt_SE/2)+5;
im_recon_show = cat(4,sos(im_recon_GE(:,:,Echo_forPhase_GE,T_to_show_GE),3),sos(im_recon_SE(:,:,Echo_forPhase_SE,T_to_show_SE),3));
% figure; imshow3(im_recon_show(end:-1:1,:,:,3:3:end),[0 20],[2,3]);
figure; imshow3(im_recon_show(end:-1:1,:,:,:),[0 20],[3,3]);
%% save data
im_recon = single(im_recon);
res_a = single(squeeze(res_a));
save('Recon_allimages','res_a','im_recon','-v7.3');
