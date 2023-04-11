h = inline('[x(1,:)*x(2,:); x(1,:)/x(2,:)]');
x = ndist([1;1], eye(2));
ytt1 = tt1eval(x, h);
ytt2 = tt2eval(x, h);
yut = uteval(x,h);
ymc = mceval(x, h);
plot2(ytt1, ytt2, yut, ymc, 'legend', '', 'col', 'bgrk')