clear; clc;close all

N    = 2;
ne   = 800;
nt   = 40000;
isnap= 20

f = fopen("OUTPUT/snapshots_sigma_analytical.bin","r");
u = fread(f,"float64");
u = reshape(u,[],nt/isnap);


figure()
for i=1:nt/isnap
    plot(u(:,i));
    pause(.02)
end