sm=exsensor('gps2d', 2,1);
sm.x0 = [1,1];
sm.pe=0.01*eye(4);
% plot(sm)
sm
y=simulate(sm, 1)
lh2(sm, y)
[xls, sls]=ls(sm,y)
[xwls, swls]=wls(sm,y)
plot(sm, sls, swls)
xplot2(xls,xwls,'conf',90)