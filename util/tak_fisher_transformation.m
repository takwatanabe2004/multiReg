function z = tak_fisher_transformation(r)
%| fisher transformation of pearon's rcorr
z = 0.5*log((1+r)./(1-r));

%| 03/17/2013
%|  I experienced some case where some of the rcorr values were +1/-1...which
%|      gives +/- inf output, which causes problems since many function gives an 
%|      error when any single element of an array containts inf....
%|      my current hacky solution: set the +/-inf value to zero (since i don't 
%|      want to use this anamolous feature for training/classification)
if sum(isinf(z(:)))~=0
    numInf = sum(isinf(z(:)));
    warning(['There were ',num2str(numInf), ' inf values...these are set to 0'])
    z(isinf(z))=0;
end