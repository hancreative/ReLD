function [ I_noisy ] = AddNoise( I, NoiseModel, sigma, imSize , percent)


if strcmp(NoiseModel, 'Gaussian')
    
    randn('seed', 0);
    noise = randn(size(I))*sigma;
    I_noisy = I + noise;
end

if strcmp(NoiseModel, 'SaltPepper')
        for t = 1:size(I,2)
        
        temp = uint8(reshape(I(:,t), imSize));
        temp=imnoise(temp,'salt & pepper',percent);
        I_noisy(:,t)=double(reshape(temp,[size(I,1),1]));
        end
end

if strcmp(NoiseModel, 'Mixed')
    
    for t = 1:size(I,2)
        
        temp = uint8(reshape(I(:,t), imSize));
        temp=imnoise(temp,'salt & pepper',percent);
        I_noisy(:,t)=double(reshape(temp,[size(I,1),1]));
    end

    randn('seed', 0);
    noise = randn(size(I))*sigma;
    I_noisy = I_noisy + noise;
end



end

