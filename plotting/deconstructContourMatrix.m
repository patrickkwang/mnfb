function [samples,level] = deconstructContourMatrix(C)
i = 1;
samplesUsed = 0;
while samplesUsed<size(C,2)
	level(i) = C(1,samplesUsed+1);
	nPoints(i) = C(2,samplesUsed+1);
	samples{i} = C(:,samplesUsed+1+(1:nPoints(i)));
	i = i+1;
	samplesUsed = sum(nPoints)+length(nPoints);
end