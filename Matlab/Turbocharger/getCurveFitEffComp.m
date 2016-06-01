function [p, resnorm] = getCurveFitEffComp(SAECompMap, polyOrder)
%%
xData = [SAECompMap.PR SAECompMap.m_dot];
yData = SAECompMap.eff;

switch polyOrder 
    case 3
        p0 = zeros(10,1);
        [p, resnorm] = lsqcurvefit(@myFun3,p0,xData,yData); 
    case 4
        p0 = zeros(15,1);
        [p, resnorm] = lsqcurvefit(@myFun4,p0,xData,yData);
    case 5
        p0 = zeros(21,1);
        [p, resnorm] = lsqcurvefit(@myFun5,p0,xData,yData);        
end;
        


x = 0:max(SAECompMap.m_dot)/100:max(SAECompMap.m_dot);
y = 1:(max(SAECompMap.PR)-1)/100:max(SAECompMap.PR);
[x, y] = meshgrid(x,y);
eff_fit = getEff(p,y,x,polyOrder);
figure
contour(x,y,eff_fit,[0.6;0.65;0.7;0.72;0.74;0.76;0.78;0.80;0.82;0.84;0.86]); 
hold on
plot3(xData(:,2),xData(:,1),yData,'o');
figure
surf(x,y,eff_fit); 
hold on
plot3(xData(:,2),xData(:,1),yData,'o');
zlim([0 1]);


end
function F = myFun3(p,xData)
x1 = xData(:,1);
x2 = xData(:,2);
F = p(1) + x2.*(p(2) + x2.*(p( 3) + x2.*p( 4))) + ...
           x1.*(p(5) + x2.*(p( 6) + x2.*p( 7)) + ...
                       x1.*(p( 8) + x2.*p( 9) + ...
                                    x1.*p(10))); 

end
function F = myFun4(p,xData)
x1 = xData(:,1);
x2 = xData(:,2);
F = p(1) + x2.*(p(2) + x2.*(p( 3) + x2.*(p( 4) + x2.*p( 5)))) + ...
           x1.*(p(6) + x2.*(p( 7) + x2.*(p( 8) + x2.*p(9))) + ...
                       x1.*(p(10) + x2.*(p(11) + x2.*p(12)) + ...
                                    x1.*(p(13) + x2.*p(14) + ...
                                                 x1.*p(15))));

end
function F = myFun5(p,xData)
x1 = xData(:,1);
x2 = xData(:,2);
F = p(1) + x2.*(p(2) + x2.*(p( 3) + x2.*(p( 4) + x2.*(p( 5) + x2.*p( 6))))) + ...
           x1.*(p(7) + x2.*(p( 8) + x2.*(p( 9) + x2.*(p(10) + x2.*p(11)))) + ...
                       x1.*(p(12) + x2.*(p(13) + x2.*(p(14) + x2.*p(15))) + ...
                                    x1.*(p(16) + x2.*(p(17) + x2.*p(18)) + ...
                                                 x1.*(p(19) + x2.*p(20) + ...
                                                              x1.*p(21)))));

end
function F = getEff(p,x1,x2,polyOrder)
switch polyOrder
    case 3
        F = p(1) + x2.*(p(2) + x2.*(p( 3) + x2.*p( 4))) + ...
                   x1.*(p(5) + x2.*(p( 6) + x2.*p( 7)) + ...
                               x1.*(p( 8) + x2.*p( 9) + ...
                                            x1.*p(10)));         
    case 4
        F = p(1) + x2.*(p(2) + x2.*(p( 3) + x2.*(p( 4) + x2.*p( 5)))) + ...
                   x1.*(p(6) + x2.*(p( 7) + x2.*(p( 8) + x2.*p(9))) + ...
                               x1.*(p(10) + x2.*(p(11) + x2.*p(12)) + ...
                                            x1.*(p(13) + x2.*p(14) + ...
                                                         x1.*p(15))));        
    case 5
        F = p(1) + x2.*(p(2) + x2.*(p( 3) + x2.*(p( 4) + x2.*(p( 5) + x2.*p( 6))))) + ...
                   x1.*(p(7) + x2.*(p( 8) + x2.*(p( 9) + x2.*(p(10) + x2.*p(11)))) + ...
                               x1.*(p(12) + x2.*(p(13) + x2.*(p(14) + x2.*p(15))) + ...
                                            x1.*(p(16) + x2.*(p(17) + x2.*p(18)) + ...
                                                         x1.*(p(19) + x2.*p(20) + ...
                                                                      x1.*p(21)))));
end;

end
%%
