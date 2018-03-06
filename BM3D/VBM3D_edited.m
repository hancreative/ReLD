function [BM3D_out ] = VBM3D_new2( sigma, TotalFrames, imSize, Data )


z = zeros([imSize(1), imSize(2), TotalFrames], 'single');
OffSet = min(0, min(Data(:)));

Data = Data + abs(OffSet)*ones(size(Data));
MaxVal = max(max(Data(:)),255);

for cf = 1:TotalFrames
        z(:,:,cf) = single(reshape(Data(:,cf),imSize)) / MaxVal; % 1/255 = 0.0039216
end
bm3dProfile         = 'np';

transform_2D_HT_name     = 'bior1.5'; %% transform used for the HT filt. of size N1 x N1
transform_2D_Wiener_name = 'dct';     %% transform used for the Wiener filt. of size N1_wiener x N1_wiener
transform_3rd_dim_name   = 'haar'; %% tranform used in the 3-rd dim, the same for HT and Wiener filt.

%%%% Step 1: Hard-thresholding (HT) parameters:
denoiseFrames       = min(9, TotalFrames); % number of frames in the temporalwindow (should not exceed the total number of frames 'NumberOfFrames')
N1                  = 8;  %% N1 x N1 is the block size used for the hard-thresholding (HT) filtering
Nstep               = 6;  %% sliding step to process every next refernece block
N2                  = 8;  %% maximum number of similar blocks (maximum size of the 3rd dimension of the 3D groups)
Ns                  = 7;  %% length of the side of the search neighborhood for full-search block-matching (BM)
Npr                 = 5;  %% length of the side of the motion-adaptive search neighborhood, use din the predictive-search BM
tau_match           = 3000; %% threshold for the block distance (d-distance)
lambda_thr3D        = 2.7; %% threshold parameter for the hard-thresholding in 3D DFT domain
dsub                = 7;  %% a small value subtracted from the distnce of blocks with the same spatial coordinate as the reference one 
Nb                  = 2;  %% number of blocks to follow in each next frame, used in the predictive-search BM
beta                = 2.0; %% the beta parameter of the 2D Kaiser window used in the reconstruction

%%%% Step 2: Wiener filtering parameters:
denoiseFramesW      = min(9, TotalFrames);
N1_wiener           = 7;
Nstep_wiener        = 4;
N2_wiener           = 8;
Ns_wiener           = 7;
Npr_wiener          = 5;
tau_match_wiener    = 1500;
beta_wiener         = 2.0;
dsub_wiener         = 3;
Nb_wiener           = 2;

%%%% Block-matching parameters:
stepFS              = 1; %% step that forces to switch to full-search BM, "1" implies always full-search
smallLN             = 3; %% if stepFS > 1, then this specifies the size of the small local search neighb.
stepFSW             = 1;
smallLNW            = 3;
thrToIncStep        = 8;  %% used in the HT filtering to increase the sliding step in uniform regions


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Parameters for the Low Complexity Profile.
%%%%
if strcmp(bm3dProfile, 'lc') == 1,
    lambda_thr3D = 2.8;
    smallLN   = 2;
    smallLNW  = 2;
    denoiseFrames  = min(5, NumberOfFrames);
    denoiseFramesW = min(5, NumberOfFrames);
    N2_wiener = 4;
    N2 = 4;
    Ns = 3;
    Ns_wiener = 3;
    NB = 1;
    Nb_wiener = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Parameters for the High Profile.
%%%%
if strcmp(bm3dProfile, 'hi') == 1,
    Nstep        = 3;
    Nstep_wiener = 3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Parameters for the "Very Noisy" Profile.
%%%%
if sigma > 30,
    N1 = 8;
    N1_wiener = 8;
    Nstep = 6;
    tau_match    = 4500;
    tau_match_wiener    = 3000;
end

decLevel           = 0;    %% dec. levels of the dyadic wavelet 2D transform for blocks (0 means full decomposition, higher values decrease the dec. number)
decLevel3          = 0;    %% dec. level for the wavelet transform in the 3rd dimension

[Tfor, Tinv]       = getTransfMatrix(N1, transform_2D_HT_name, decLevel); %% get (normalized) forward and inverse transform matrices
[TforW, TinvW]     = getTransfMatrix(N1_wiener, transform_2D_Wiener_name); %% get (normalized) forward and inverse transform matrices
thr_mask           = ones(N1); %% N1xN1 mask of threshold scaling coeff. --- by default there is no scaling, however the use of different thresholds for different wavelet decompoistion subbands can be done with this matrix

if (strcmp(transform_3rd_dim_name, 'haar') == 1 || strcmp(transform_3rd_dim_name(end-2:end), '1.1') == 1),
    %%% Fast internal transform is used, no need to generate transform
    %%% matrices.
    hadper_trans_single_den         = {};
    inverse_hadper_trans_single_den = {};
else
    %%% Create transform matrices. The transforms are later computed by
    %%% matrix multiplication with them
    for hh = [1 2 4 8 16 32];
        [Tfor3rd, Tinv3rd]   = getTransfMatrix(hh, transform_3rd_dim_name, decLevel3);
        hadper_trans_single_den{hh}         = single(Tfor3rd);
        inverse_hadper_trans_single_den{hh} = single(Tinv3rd');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 2D Kaiser windows that scale the reconstructed blocks
%%%%
if beta_wiener==2 & beta==2 & N1_wiener==7 & N1==8 % hardcode the window function so that the signal processing toolbox is not needed by default
    Wwin2D = [ 0.1924    0.2989    0.3846    0.4325    0.4325    0.3846    0.2989    0.1924;
        0.2989    0.4642    0.5974    0.6717    0.6717    0.5974    0.4642    0.2989;
        0.3846    0.5974    0.7688    0.8644    0.8644    0.7688    0.5974    0.3846;
        0.4325    0.6717    0.8644    0.9718    0.9718    0.8644    0.6717    0.4325;
        0.4325    0.6717    0.8644    0.9718    0.9718    0.8644    0.6717    0.4325;
        0.3846    0.5974    0.7688    0.8644    0.8644    0.7688    0.5974    0.3846;
        0.2989    0.4642    0.5974    0.6717    0.6717    0.5974    0.4642    0.2989;
        0.1924    0.2989    0.3846    0.4325    0.4325    0.3846    0.2989    0.1924 ];
    Wwin2D_wiener = [ 0.1924    0.3151    0.4055    0.4387    0.4055    0.3151    0.1924;
        0.3151    0.5161    0.6640    0.7184    0.6640    0.5161    0.3151;
        0.4055    0.6640    0.8544    0.9243    0.8544    0.6640    0.4055;
        0.4387    0.7184    0.9243    1.0000    0.9243    0.7184    0.4387;
        0.4055    0.6640    0.8544    0.9243    0.8544    0.6640    0.4055;
        0.3151    0.5161    0.6640    0.7184    0.6640    0.5161    0.3151;
        0.1924    0.3151    0.4055    0.4387    0.4055    0.3151    0.1924 ];
else
    Wwin2D           = kaiser(N1, beta) * kaiser(N1, beta)'; % Kaiser window used in the aggregation of the HT part
    Wwin2D_wiener    = kaiser(N1_wiener, beta_wiener) * kaiser(N1_wiener, beta_wiener)'; % Kaiser window used in the aggregation of the Wiener filt. part
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Read an image, generate noise and add it to the image
%%%%

l2normLumChrom = ones(TotalFrames,1); %%% NumberOfFrames == nSl !

% if dump_information == 1,
%     fprintf('Video: %s (%dx%dx%d), sigma: %.1f\n', Xnoisy_name, videoHeight, videoWidth, NumberOfFrames, sigma);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initial estimate by hard-thresholding filtering

y_hat = bm3d_thr_video(z, hadper_trans_single_den, Nstep, N1, N2, 0,...
    lambda_thr3D, tau_match*N1*N1/(MaxVal*MaxVal), (Ns-1)/2, sigma/MaxVal, thrToIncStep, single(Tfor), single(Tinv)', inverse_hadper_trans_single_den, single(thr_mask), 'unused arg', dsub*dsub/MaxVal, l2normLumChrom, Wwin2D, (Npr-1)/2, stepFS, denoiseFrames, Nb );



y_hat_wi = bm3d_wiener_video(z, y_hat, hadper_trans_single_den, Nstep_wiener, N1_wiener, N2_wiener, ...
    'unused_arg', tau_match_wiener*N1_wiener*N1_wiener/(MaxVal*MaxVal), (Ns_wiener-1)/2, sigma/MaxVal, 'unused arg', single(TforW), single(TinvW)', inverse_hadper_trans_single_den, 'unused arg', dsub_wiener*dsub_wiener/MaxVal, l2normLumChrom, Wwin2D_wiener, (Npr_wiener-1)/2, stepFSW, denoiseFramesW, Nb_wiener );

% In case the input noisy video is clipped in [0,1], then apply declipping  
% if isCharacterName
%     if exist('Xorig', 'var') == 1
%         if ~strcmp(Xorig, Xnoisy_name)
             [y_hat_wi] = ClipComp16b(sigma/MaxVal, y_hat_wi);
%         end
%     else
%         [y_hat_wi] = ClipComp16b(sigma/255, y_hat_wi);
%     end
% end



y_hat_wi = y_hat_wi*MaxVal;
y_hat_wi = y_hat_wi-abs(OffSet)*ones(size(y_hat_wi));

for i = 1:size(y_hat_wi,3)
 temp = y_hat_wi(:,:,i);
 temp = double(temp);
 
 BM3D_out(:,i) = reshape(temp,[imSize(1)*imSize(2),1] );
end


end

