de = [1e-7 5e-8 1e-8 5e-9 1e-9 5e-10 1e-10 1e-11];
for i = 5:8
    disp(de(i))
   p.D_e = de(i);
   run dfn_scott_PSO_June_2016
end