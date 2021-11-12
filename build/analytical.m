clear; clc; close all;
nt   = 60000;
ne   = 800;
ircv = 100;
nrcv = 8; 
N    = 5;
h    = 5.;
esrc = 400;
gsrc = 1;
dt   = .000015;
f0   = 6.;

t0   = dt/f0;

f   = fopen("OUTPUT/source.bin","r");
src = fread(f,"float64");
fclose(f);

f   = fopen("testing_vs","r");
v   = fread(f,"float64");
fclose(f);

c   = max(v);


xi= [-1.0, -0.765055323, -0.285231516,   0.28523151, 0.765055323, 1.0 ];     
xg   = zeros(ne,N+1);

for i=1:ne
    for j=1:N+1
            xg(i,j) = h*(i-1) +  h * ((xi(j) + 1) / 2);
    end
end 



% Green's function


greensfn = zeros(nt,nrcv);      
xs = xg(esrc,gsrc);

for ir=1:nrcv
xr  = xg(ir*ircv,ceil((N+1)/2));
    for i=1:nt
        if (((i*dt-t0)-abs(xr-xs)/c)>=0)
            greensfn(i,ir)=1./(2*c);
        else
            greensfn(i,ir)=0;
        end
    end
end

greens = zeros(nt,nrcv);      


for ir=1:nrcv
    greens(:,ir) = conv(greensfn(:,ir),src,"same");
end


plot(conv(greensfn(:,1),src))