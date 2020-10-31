clear all, close all, clc

energy1 = importdata('Energyequilibrium.dat');
energy2 = importdata('Energyequilibrium2.dat');

hold on

plot(energy1(:,1), energy1(:,4), 'Linewidth',1.5);     %%TE
 plot(energy2(:,1), energy2(:,4), 'Linewidth',1.5);     %%TE

%  legend('PE', 'KE', 'TE');



