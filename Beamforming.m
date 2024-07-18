%均匀线阵的阵元间距取半波长
clc;
clear;
close all;

[theta, G_5] = Calc_G(5);
[~, G_10] = Calc_G(10);
[~, G_15] = Calc_G(15);
[~, G_20] = Calc_G(20);

figure('Color','white','Name','frequency-wavenumber response','NumberTitle','off','Position', [500,500,800,480])
subplot(2,2,1); polardb(theta,G_5,-40,'b-'); title('N=5')
subplot(2,2,2); polardb(theta,G_10,-40,'b-');title('N=10')
subplot(2,2,3); polardb(theta,G_15,-40,'b-');title('N=15')
subplot(2,2,4); polardb(theta,G_20,-40,'b-');title('N=20')


function [theta, G_dB] = Calc_G(N)
    n = (0 : N - 1);
    lambda = 1;%波长取单位长度
    d = lambda / 2;
    theta = (0 : 0.001 : 2 * pi);
    V_theta = exp(-1i * (n - (N - 1) / 2)' * 2 * pi / lambda * cos(theta) * d);
    G_dB = 20 * log10(abs(ones(1, N) * V_theta / N));
end

