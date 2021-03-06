function [mDotCorr, PR, effIs] = ...
    convCorrectedComp(mDot,pIn,pOut,TIn,TOut,FIn)
mDotCorr = mDot*1e5./pIn.*sqrt(TIn/298);
PR = pOut./pIn;
[~,hIn,~,~,~,~,~,~,~,~,~,~,~,CpIn,~,KIn] = ...
    GetThdynCombGasZachV1(pIn,TIn,FIn,0.0683);
[~,hOut,~,~,~,~,~,~,~,~,~,~,~,CpOut,~,~] = ...
    GetThdynCombGasZachV1(pOut,TOut,FIn,0.0683);
TIs = TIn.*PR.^((KIn-1)/KIn);
[~,hOutIs,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
    GetThdynCombGasZachV1(pOut,TIs,FIn,0.0683);

effIs = (hOutIs - hIn)./(hOut - hIn);