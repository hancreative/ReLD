function sigma = StdEst( data )
numTrain = min(100, size(data,2));
SigmaTrain = data(:,1:numTrain);
mu = mean(SigmaTrain,2);
MTrain=SigmaTrain-repmat(mu,1,numTrain);
sigma = sqrt(var(MTrain(:)));

end

