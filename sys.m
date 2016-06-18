function dy = sys(t, y)
alpha = 4;
gamma = 1.1;
dy = [ y(1).*(alpha - y(1) - y(2)./(1 + y(1))); y(2).*(-1 + gamma*y(1)./(1 + y(1)))];
