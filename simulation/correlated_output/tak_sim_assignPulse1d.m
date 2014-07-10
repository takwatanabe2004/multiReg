function w = tak_sim_assignPulse1d(p, nPatches,patchLenBase)
% w = tak_sim_assignPulse1d(p, nPatches)
%=========================================================================%
% - Assign "nPatches" of "nonoverlapping" smooth pulses on a 1d signal 
% (patchLenBase = "base-length" of each patches)
%=========================================================================%
% (07/07/2014)
%%
w=zeros(p,1);
if nargin==2
    patchLenBase= round(p/100);
end
for ipatch = 1:nPatches
    %=====================================================================%
    % create arbitrary patch for each outputs
    %=====================================================================%   
    patchLen = patchLenBase * (3+randsample(3,1)); % patchLen \in K*{4,5,6}
    
    % frequency (if numerator here = 2, it means full sinusoidal cycle)       
    freq = randsample(6,1)/(1*patchLen);
        
    % magnitude of sinusoidal pulses
    mag = tak_sample_signed_unif([5,12],1);
        
%     mag=abs(mag);
    %=====================================================================%
    % ensure patches don't overlap
    % - assign patch at interval region (1:p/npatch) + offset
    % - so length of each patch is less than this interval width
    % (interval width = p/npatch)
    %=====================================================================%
    intervalWidth = round(p/nPatches);
    offset = intervalWidth*(ipatch-1);
    patchStart = offset+randsample(intervalWidth-patchLen-1,1); %
    
    %         W(patchStart+1:patchStart+patchLen,iq) = tak_sample_signed_unif([5,10],1);
    tmp = zeros(p,1);
    tmp(patchStart+1:patchStart+patchLen) = mag*sin(freq*pi*(1:patchLen));
%     keyboard
    w = w + tmp;
end
