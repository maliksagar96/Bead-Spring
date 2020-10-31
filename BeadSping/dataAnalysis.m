clear all, close all, clc

energy = importdata('Energy.dat');

hold on

% plot(energy(:,1), energy(:,2), 'Linewidth',1.5);     %%PE
% plot(energy(:,1), energy(:,3), 'Linewidth',1.5);     %%KE
plot(energy(:,1), energy(:,4), 'Linewidth',1.5);     %%TE

%  legend('PE', 'KE', 'TE');



