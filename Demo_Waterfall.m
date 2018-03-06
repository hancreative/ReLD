clear all; clc;
addpath ReProCS;
addpath ReProCS/Yall1;
addpath inexact_alm_rpca;
addpath inexact_alm_rpca/PROPACK;
addpath Data;
addpath BM3D;


load Waterfall_Small.mat;

sigma = 70;

NoiseModel = 'Gaussian';
I_noisy = AddNoise(I, NoiseModel, sigma, imSize);

%denoise using ReLD
[Shat_repro, Lhat_repro] = ReProCS(I_noisy, imSize);
Ldenoised_repro = VBM3D_edited(StdEst(Lhat_repro), size(Lhat_repro,2), imSize, Lhat_repro);
Sdenoised_repro = VBM3D_edited(StdEst(Shat_repro), size(Shat_repro,2), imSize, Shat_repro);
Idenoised_repro = Ldenoised_repro + Sdenoised_repro;


%compare with denoising using VBM3D;
VBM3D_denoised = VBM3D_edited(StdEst(I_noisy), size(I_noisy,2), imSize, I_noisy);

PSNR_ReLD = 20*log10(255 * sqrt(numel(I_noisy)) / norm(Idenoised_repro(:)-I(:)));
PSNR_VBM3D = 20*log10(255 * sqrt(numel(I_noisy)) / norm(VBM3D_denoised(:)-I(:)));

DisplayVideo(I, I_noisy, Idenoised_repro, VBM3D_denoised, imSize, 'Waterfall_result');