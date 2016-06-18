%% Неподвижные точки и устойчивость
cla reset
u0 = 0.00001;
r = 2.5; % r = 2;
t = 20;
ut = zeros(1, t);
ut1 = zeros(1, t);
u_prev = u0;
repel = 0;

attract = (1 + 2*r*r - sqrt(1 + 4 * r * r))/(2* r* r);
for i = 1 : t
    ut(i) = u_prev;   
    %ut1(i) = sqrt(u_prev)*exp(r*(1 - u_prev*u_prev));
    ut1(i) = r *sqrt(u_prev)*(1 - u_prev);
    u_prev = ut1(i);
end
%plot(ut, ut1,'b')
f=@(u) r.*sqrt(u).*(1-u);
net=linspace(0, 1, 500);
f_net=f(net);
plot(net, f_net, 'b');
xlabel('u');
ylabel('f(u)');
hold on
plot(ut, ut, 'g')
plot(repel, repel, '*y');
plot(attract, attract, '*r');

point = attract - 0.2;
x0 = point;
y0 = 0;
while (abs(point - attract) > 0.01)
    next =  r *sqrt(point)*(1 - point);
    x1 = point;
    y1 = next;
    quiver(x0, y0, x1 - x0, y1 - y0, 0, 'k');
    y0 = next;
    x1 = next; 
    quiver(x0, y0, x1 - x0, y1 - y0, 0, 'k');
    x0 = x1;
    point = next;
end
point = attract + 0.2;
x0 = point;
y0 = 0;
while (point ~= attract)
    next =  r *sqrt(point)*(1 - point);
    x1 = point;
    y1 = next;
    quiver(x0, y0, x1 - x0, y1 - y0, 0, 'c');
    y0 = next;
    x1 = next; 
    quiver(x0, y0, x1 - x0, y1 - y0, 0, 'c');
    x0 = x1;
    point = next;
end
point = repel + 0.06;
x0 = point;
y0 = 0;
i = 0;
while (i < 3)
    i = i + 1;
    next =  r *sqrt(point)*(1 - point);
    x1 = point;
    y1 = next;
    quiver(x0, y0, x1 - x0, y1 - y0, 0, 'm');
    y0 = next;
    x1 = next; 
    quiver(x0, y0, x1 - x0, y1 - y0, 0, 'm');
    x0 = x1;
    point = next;
    i
end

print('nepodv','-dpng')
ut;
ut1;
%% Бифуркационная диаграмма
r = 1 : 0.02 : 6;
hold on
n1 = 200;
n2 = 100;
u0 = 0.01;
for i = 1 : length(r)
    u_prev = u0;
    for j = 1 : n1
  
        ut1 = r(i) *sqrt(u_prev)*(1 - u_prev);
        u_prev = ut1;
    end
    for l = j : n1 + n2
  
        ut1 = r(i) *sqrt(u_prev)*(1 - u_prev);
        plot(r(i), ut1, 'b.');
        u_prev = ut1;
    end     
end
xlabel('r');
ylabel('u_t,');
print('bifurk','-dpng');

%% Цикл длины 3
clear
r = 2.59;
f = @(u) r.*sqrt(u).*(1 - u);
g = @(u) f(f(f(u)));
syms x;
MyFun(x)=sym(g(x));
p = diff(MyFun);
u = 0.2 : 0.007 : 1 - 0.01;
for i = 1 : length(u)
    w(i) = p(u(i));
end
plot(u, g(u) - u, 'g');

hold on
plot(u, w - 1, 'r');
xlabel('u');
ylabel('function');
legend('f(f(f(u)))-u', 'df(f(f(u))) - 1');
print('r3','-dpng');
%% Цикл график
cla reset
u0 = 0.00001;
r = 2.5815; % r = 2;
t = 20;
ut = zeros(1, t);
ut1 = zeros(1, t);
u_prev = u0;
repel = 0;

attract = (1 + 2*r*r - sqrt(1 + 4 * r * r))/(2* r* r);
for i = 1 : t
    ut(i) = u_prev;   
    %ut1(i) = sqrt(u_prev)*exp(r*(1 - u_prev*u_prev));
    ut1(i) = r *sqrt(u_prev)*(1 - u_prev);
    u_prev = ut1(i);
end
%plot(ut, ut1,'b')
f=@(u) r.*sqrt(u).*(1-u);
net=linspace(0, 1, 500);
f_net=f(net);
plot(net, f_net, 'b');
xlabel('u');
ylabel('f(u)');
hold on
plot(ut, ut, 'g')
plot(repel, repel, '*y');
plot(attract, attract, '*r');

point = 0.33;
x0 = point;
y0 = 0;
iter = 5;
i = 0;
while (i < iter)
    i = i + 1;
    next =  r *sqrt(point)*(1 - point);
    x1 = point;
    y1 = next;
    quiver(x0, y0, x1 - x0, y1 - y0, 0, 'k');
    y0 = next;
    x1 = next; 
    quiver(x0, y0, x1 - x0, y1 - y0, 0, 'k');
    x0 = x1;
    point = next;
end


print('cycle','-dpng')

%% Показатель Ляпунова
cla reset
r = 0.01 : 0.001 : 4;
hold on
n = 10000;
u0 = 0.5;
u_prev = u0;
f = @(u)sqrt(u)*(1 - u);
p = @(u, r)r*(1/sqrt(u) - 3*sqrt(u))/2;
for i = 1 : length(r)
    lyap = 0;
    u_prev = u0;
    for l = 1 : n
        ut = u_prev;   
        ut1 = r(i) *f (u_prev);
  
        y = abs(p(r(i),ut));
        lyap = lyap + log(y);
        u_prev = ut1;
    end 
    plot(r(i), lyap/n, 'b.');
   
end
plot(r, 0.*r,'c')
xlabel('r');
ylabel('h');
print('lyapunov','-dpng');
%% Нейман - Сакер
r = sqrt(2) + 0.01
iter = 10000;
n = iter;
u = zeros(1, n);
v = zeros(1, n);
u(1) = 0.2;
u(2) = 0.2;
v(1) = 0.3;
v(2) = 0.3;
i = 2;
while i < iter
    u(i + 1) = r * sqrt(u(i)).*(1 - u(i - 1));
    v(i) = u(i - 1);
    i = i + 1;
end
plot(v, u, 'b.');
xlabel('v');
ylabel('u');
print('saker2+','-dpng');
%% Фазовый портрет
 set(0,'DefaultTextInterpreter', 'latex')
gamma = 1.05: 0.01 : 1.5;
a1 = @(g)1./(g - 1);
a2 = @(g)1 + 2./(g - 1) + 2*g.*(g - 1) - 2*g.*sqrt((g - 1).^2 + 1);
a3 = @(g)1 + 2./(g - 1);
a4 = @(g)1 + 2./(g - 1) + 2*g.*(g - 1) + 2*g.*sqrt((g - 1).^2 + 1);
plot(gamma, a1(gamma), 'r');
hold on
plot(gamma, a2(gamma), 'g');
plot(gamma, a3(gamma), 'b');
plot(gamma, a4(gamma), 'm');
xlabel('\gamma');
ylabel('\alpha');
legend('1', '2', '3', '4', '5', '6');
text(1.1, 5, 'A');
text(1.1, 12, 'B');
text(1.2, 9.5, 'C');
text(1.3, 10, 'D');
text(1.3, 15, 'E');
print('faz_p', '-dpng')
%% Фазовый портрет системы (x, y)
x0 =0: 0.35 : 10;
y0 = 0 : 0.35: 10;
for i = 1 : numel(x0)
    for j = 1 : numel(y0)
        [t y] = ode45(@sys, [0 1], [x0(i) y0(j)]);
    plot(y(:,1), y(:,2), 'b')
    hold on
    end
end
alpha = 4;
gamma = 1.1;
Mx = 1./(-1 + gamma);
My = (alpha - 1./(-1 + gamma)).*(1 + 1./(-1 + gamma));
if ((Mx >= 0)&&(My >=0))
    plot(Mx, My, 'r*');
    text(Mx, My, 'M');
end
plot(0, 0, 'r*');
    text(0, 0, 'O');
plot(alpha, 0, 'r*');
    text(alpha, 0, 'K');
xlabel('x');
ylabel('y');
print('fp1','-dpng');
