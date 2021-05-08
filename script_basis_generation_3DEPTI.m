% Example script of the basis generation for 3D-EPTI subspace reconstruction
% Related to the following works
% 1) "Ultrafast high-resolution multi-parametric brain MRI using 3D Echo Planar Time-resolved Imaging",Fuyixue Wang and Zijing Dong et al., 2021.
% 2) "Echo planar time-resolved imaging (EPTI)", Fuyixue Wang et al.,MRM,2019
% 3) "EPTI with Subspace Reconstruction and Optimized Spatiotemporal Encoding", Zijing Dong et al.2020
% Fuyixue Wang <fwang18@mgh.harvard.edu>, Zijing Dong <zijingd@mit.edu> Feb/2021

clear;
close all;
addpath(genpath('Funcs'));
%% simulation parameters
load('acq_params.mat')

T2s=[10:1:100,102:2:200,204:10:400,420:20:500]';
T2 = T2s;
T1=[400:200:2000,2400:500:5000]';

N1=128; % maximum number of T2* value used
N2=128; % maximum number of T2 value used

B1_factor = 0.75:0.05:1.25;
% B1_factor=0.8;
T2s=T2s/1000;
T2=T2/1000;
T1=T1/1000;

K=12;
flag_showfig=1;

[U, X] = gen_basis_withB1(N1, N2, T2s, T2, T1, B1_factor, param);
Phi = U(:,1:K);
Z = Phi*Phi'*X;
err = norm(X(:) - Z(:)) / norm(X(:));
fprintf('Relative norm of error: %.6f\n', err);
%%
TEs_all=param.Time_sampling(1:param.N_signal);
if flag_showfig==1
    figure;
    plot(TEs_all*1000, real(X(:,1:1000:end)), 'linewidth', 2); hold on;
    xlabel('Virtual echo time (ms)');
    ylabel('Signal value');
    ftitle('Signal evolutions for distribution of T2 T2* and T1 values', 24)
    faxis;

    figure;
    subplot(1,2,1); plot(real(Phi), 'linewidth', 3);
    ftitle('Subspace curves Real Part', 24)
    faxis;
    subplot(1,2,2); plot(imag(Phi), 'linewidth', 3);
    ftitle('Subspace curves Imag Part', 24)
    faxis;
    %% Project the signal evolutions onto the subspace
    figure;
    plot(TEs_all*1000, real(Z(:,1:1000:end)), 'linewidth', 2); hold on;
    xlabel('Virtual echo time (ms)');
    ylabel('Signal value');
    ftitle('Projected signal evolutions', 24)
    faxis;
end
save(['Generated_bases.mat'],'U','-v7.3');
disp('Done saving data!');
