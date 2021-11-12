clear; clc; close all;

nt   = 60000;
ne   = 800;
ircv = 100;

f = fopen("OUTPUT/seism_sigma.bin","r");
u = fread(f,"float64");
fclose(f);

u = reshape(u,nt,[]);