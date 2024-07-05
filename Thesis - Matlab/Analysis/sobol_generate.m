clear;clc;close all;
boundaries = [-10,10;-10,10];
dim = 2;
npoint = 100;
sob_point = sobolgen(npoint,dim,boundaries)

% Generate titik-titik acak dengan distribusi seragam
x = unifrnd(boundaries(1,1), boundaries(1,2), 1, npoint);
y = unifrnd(boundaries(2,1), boundaries(2,2), 1, npoint);

% Plot
figure;
fz = 28;
xlim(boundaries(1,:));
ylim(boundaries(2,:));
subplot(1,2,1);
scatter(x, y, 'filled','MarkerFaceColor','#224141');
title('(a) Distribusi Bilangan Acak','FontSize',fz);
xlabel('x');
ylabel('y');
grid on;
% axis equal;
pbaspect([diff(xlim()) diff(ylim()) 1]);

subplot(1,2,2);
scatter(sob_point(:,1), sob_point(:,2), 'filled','MarkerFaceColor','#B185B4');
title('(b) Distribusi Bilangan Sobol','FontSize',fz);
xlabel('x');
ylabel('y');
grid on;
% axis equal;
pbaspect([diff(xlim()) diff(ylim()) 1]);

set(gcf, 'WindowState', 'maximized');