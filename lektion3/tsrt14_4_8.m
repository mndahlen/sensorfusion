% h = inline('[x(1,:)*x(2,:); x(1,:)/x(2,:)]');
sm = sensormod('[x(1,:)*x(2,:); x(1,:)/x(2,:)]', [2,0,2,0]);
sm.x0 = [1,1];
sm.pe = ndist([0;0], [0.1, 0.05; 0.05, 0.3]);
% sm.px0 = [0,0];

y = simulate(sm, 10)
[xls, sls]=ls(sm,y)
[xwls, swls]=wls(sm,y)
plot(sm, sls, swls)
xplot2(xls,xwls,'conf',90)
