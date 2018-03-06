clear all; clc;
addpath ReProCS;
addpath ReProCS/Yall1;
addpath inexact_alm_rpca;
addpath inexact_alm_rpca/PROPACK;
addpath Data;
addpath grasta;
addpath ncrpca;
addpath ORPCA;
addpath rpca-gd;
addpath rpca-gd/private1
addpath rpca-gd/MinMaxSelection
addpath BM3D;

load Waterfall_Small.mat;
sigma = 70;
NoiseModel = 'Gaussian';

I_noisy = AddNoise(I, NoiseModel, sigma, imSize);
global rhat;
%denoise using ReLD
[Shat_repro, Lhat_repro] = ReProCS(I_noisy, imSize);
Ldenoised_repro = VBM3D_edited(StdEst(Lhat_repro), size(Lhat_repro,2), imSize, Lhat_repro);
Sdenoised_repro = VBM3D_edited(StdEst(Shat_repro), size(Shat_repro,2), imSize, Shat_repro);
Idenoised_repro = Ldenoised_repro + Sdenoised_repro;

%denoise using ORPCA

[tmp1_stoc, tmp2_stoc, S_stoc] = stoc_rpca(I_noisy, rhat);
L_stoc = I_noisy-S_stoc;
Ldenoised_stoc = VBM3D_edited(StdEst(L_stoc), size(L_stoc,2), imSize, L_stoc);
Sdenoised_stoc = VBM3D_edited(StdEst(S_stoc), size(S_stoc,2), imSize, S_stoc);
Idenoised_stoc = Sdenoised_stoc + Ldenoised_stoc;

%denoise using rpca_gd
alpha=0.1;
params.step_const = 0.5;
params.max_iter = 30;
params.tol = 2e-4;
params.incoh = 5;
[U, V] = rpca_gd(I_noisy, rhat, alpha, params);
L_gd = U*V';
S_gd = I_noisy - L_gd;
Ldenoised_gd = VBM3D_edited(StdEst(L_gd), size(L_gd,2), imSize, L_gd);
Sdenoised_gd = VBM3D_edited(StdEst(S_gd), size(S_gd,2), imSize, S_gd);
Idenoised_gd = Sdenoised_gd + Ldenoised_gd;

%denoise using non-convex rpca (AltProj)
[L_ncrpca, S_ncrpca,~,~] = ncrpca(I_noisy, rhat);
Ldenoised_ncrpca = VBM3D_edited(StdEst(L_ncrpca), size(L_ncrpca,2), imSize, L_ncrpca);
Sdenoised_ncrpca = VBM3D_edited(StdEst(S_ncrpca), size(S_ncrpca,2), imSize, S_ncrpca);
Idenoised_ncrpca = Sdenoised_ncrpca + Ldenoised_ncrpca;

%denoise using PCP
[Shat_PCP, Lhat_PCP] = PCP(I_noisy, imSize);
Ldenoised_PCP = VBM3D_edited(StdEst(Lhat_PCP), size(Lhat_PCP,2), imSize, Lhat_PCP);
Sdenoised_PCP = VBM3D_edited(StdEst(Shat_PCP), size(Shat_PCP,2), imSize, Shat_PCP);
Idenoised_PCP = Ldenoised_PCP + Sdenoised_PCP;

%denoise using GRASTA
[Shat_grasta, Lhat_grasta] = grasta(I_noisy, imSize);
Ldenoised_grasta = VBM3D_edited(StdEst(Lhat_grasta), size(Lhat_grasta,2), imSize, Lhat_grasta);
Sdenoised_grasta = VBM3D_edited(StdEst(Shat_grasta), size(Shat_grasta,2), imSize, Shat_grasta);
Idenoised_grasta = Ldenoised_grasta + Sdenoised_grasta;