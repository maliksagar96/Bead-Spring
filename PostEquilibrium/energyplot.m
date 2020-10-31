clear all, close all, clc


energy1 = importdata('Energy.dat');
energy2 = importdata('E005.dat');

hold on

plot(energy1(:,1), energy1(:,4), 'Linewidth',1.5);     %%TE
plot(energy2(:,1), energy2(:,4), 'Linewidth',1.5);     %%TE

legend('T = .025', 'T = .005');
