syms r
u = (2*r.^2 + 1 - sqrt(4*r.^2 + 1))./(2*r.^2)
u = ((sqrt(4*r.^2 + 1) - 1)./(2*r)).^2;
f = r.*(1 - 3*u)./(2*sqrt(u))
simplify(f)
r1 = solve(f == -1)
r2 = r1(2)
u1 = subs(u, r, r2)
f1 = subs(f, r, r2);
simplify(subs(f1, u, u1))