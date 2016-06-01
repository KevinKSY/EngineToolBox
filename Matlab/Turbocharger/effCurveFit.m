function p = getCurveFitEffComp(SAECompMap);
%%
xData = [SAECompMap.PR SAECompMap.m_dot];
yData = SAECompMap.eff;


[p, resnorm] = lsqcurvefit(@myFun,p0,p,xData,yData);

function F = myFun(p,xData)
x1 = xData(:,1);
x2 = xData(:,2);
F = x(1) + x2.*(p(2) + x2.*(p( 3) + x2.*(p( 4) + x2.*(p( 5) + x2.*p( 6))))) + ...
           x1.*(p(7) + x2.*(p( 8) + x2.*(p( 9) + x2.*(p(10) + x2.*p(11)))) + ...
                       x1.*(p(12) + x2.*(p(13) + x2.*(p(14) + x2.*p(15))) + ...
                                    x1.*(p(16) + x2.*(p(17) + x2.*p(18)) + ...
                                                 x1.*(p(19) + x2.*p(20)) + ...
                                                              x1.*p(21))));
end

%%
