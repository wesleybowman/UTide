addpath UTide_matlab

x = load( '-ascii', 'Honolulu/HNL2010.asc');
size(x)
t = x(:, 1);
h = x(:, 2);

t = t + datenum(1700, 1, 1, 0, 0, 0);

cnstit = 'auto'

coef = ut_solv ( t, h, [], 21, cnstit, 'OLS', 'linci');
coef_all = coef

save -v6 HNL2010.mat cnstit t h coef

n = 24 * 31
t = x(1:n, 1) + datenum(1700, 1, 1, 0, 0, 0);
h = x(1:n, 2);
coef = ut_solv(t, h, [], 21, cnstit, 'OLS', 'linci');
coef_month = coef
save -v6 HNL2010_Jan.mat cnstit t h coef

n = 24 * 7
t = x(1:n, 1) + datenum(1700, 1, 1, 0, 0, 0);
h = x(1:n, 2);
coef = ut_solv(t, h, [], 21, cnstit, 'OLS', 'linci');
save -v6 HNL2010_Jan_week1.mat cnstit t h coef

infer.infnam = {'S2  '; 'N2  '; 'O1  '};
infer.refnam = {'M2  '; 'M2  '; 'K1  '};
infer.amprat = [coef_month.A(4) / coef_month.A(2); ...
                coef_month.A(5) / coef_month.A(2); ...
                coef_month.A(3) / coef_month.A(1)];
infer.phsoff = [coef_month.g(2) - coef_month.g(4); ...
                coef_month.g(2) - coef_month.g(5); ...
                coef_month.g(1) - coef_month.g(3)];

coef = ut_solv(t, h, [], 21, cnstit, 'OLS', 'linci', 'infer', infer);
save -v6 HNL2010_Jan_week1_infer_S2.mat cnstit t h coef, infer
