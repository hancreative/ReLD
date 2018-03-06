function [ Shat, Lhat ] = PCP( Data, imSize )

lambda = 1/sqrt(max(size(Data)));
[Lhat, Shat,~] = inexact_alm_rpca(Data, lambda);


end

