function [ Shat, Lhat ] = grasta( Data, imSize )

M=Data;

global rhat;
OPTIONS.RANK                = rhat;  % the estimated low-rank
OPTIONS.rho                 = 1.8;    
OPTIONS.ITER_MAX            = 20; 
OPTIONS.ITER_MIN            = 20;    % the min iteration allowed for ADMM at the beginning

OPTIONS.USE_MEX             = 0;     % If you do not have the mex-version of Alg 2
                                     % please set Use_mex = 0.                                     

%% Initialize the rough subspace
OPTIONS.CONSTANT_STEP       = 0;   % use adaptive step-size to initialize the subspace
OPTIONS.MAX_LEVEL           = 20;
OPTIONS.MAX_MU              = 10000; % set max_mu large enough for initial subspace training
OPTIONS.MIN_MU              = 1;
FPS_ONLY                    = 1;    % 0:show the training video, 1: suppress the video demostration
TRAIN_FRAME                 = 1;    % 0¡Guse the first #training_size frames¡F 
                                    % 1: random select #training_size frames
                                    
max_cycles                  = 10;    % training cycles
training_size               = 100;   % random chose 50 frames as the training set
TRAINING_SAMPLING           = 0.3;   % Use how much information to train the first subspace.


[bgU, status, OPTS]  = bgtraining(Data, imSize, OPTIONS, max_cycles, TRAINING_SAMPLING, training_size,FPS_ONLY,TRAIN_FRAME);



%% Make video -- grasta
OPTS.MAX_ITER               = 20;
OPTIONS.CONSTANT_STEP       = 1e-2; % use the constant step-size
FPS_ONLY                    = 1;    % if you want to measure the FPS performance, please let FPS_ONLY=1
SAMPLING                    = 0.2;  % Use how much information to track the subspace.
thresh                      = 0.2;
MAX_FRAME                   = -1;   % -1 means seperating all the frames
OPTIONS.USE_MEX             = 0;
%fprintf('Seperating the whole video sequence by grasta...\n');
[video_grasta_fg,video_grasta_bg, vInfo] = bgfg_seperation_grasta( M, imSize, bgU,  SAMPLING ,status,OPTIONS, OPTS,thresh,FPS_ONLY, MAX_FRAME);

Shat = M-video_grasta_bg;
Lhat = video_grasta_bg;

end

