clear

sm = exsensor('toa', 8);
sm.th = [1,1, -1,-1, 1,-1, -1,1, 0,1, 0,-1, -1,0, 1,0];
sm.x0 = [0,0];
sm.pe = 0.1*eye(8);
y = simulate(sm, 0);
% plot(sm)

s1 = sm;
s1.x0 = [0.5,0.5];
s1.pe = 0.3*diag(rand(8,1));
% plot(s2)

xls = ls(s1, y);
xwls = wls(s1, y);
xhat = estimate(s1, y, 'thmask', zeros(s1.nn(4), 1))
% plot(s1);
% hold on;
% plot(xhat, 'conf', 90);
% hold off;

plot(sm)
hold on;
x = crlb(sm, y)
x.Px
