
clear;clc;close all;
dx = 0.5;                               % x阵元间距
dy = 0.5;                               % y间距
theta0 = 45 * (pi/180);                           % 俯仰角180
phi0 = 30 * (pi/180);                        % 方位角180
range= 10;
lamda=0.125;
alpha_x = 2*pi*sind(phi0)*cosd(theta0);      % x相位差
alpha_y = 2*pi*sind(theta0);           % y相位差        
M = 20;                                 % x阵元数
N = 20;                                 % y阵元数                      
X = (0:1:M-1)*dx;                          % x阵列
Y = (0:1:N-1)*dy;                          % y阵列
X2=kron(ones(1,N),X);
Y2=kron(Y,ones(1,M));
 
figure;
plot(X2,Y2,'.');axis equal;grid on;
title('Antenna Array');xlabel('Distance');ylabel('Distance');
 
ax = exp(1i*X*alpha_x);                 % x方向导向矢量
ay = exp(1i*Y*alpha_y);                 % y方向导向矢量
axy = kron(ax,ay);                      % 矩形阵面导向矢量
 
 
dtheta = 0.2;
dphi = 0.2;                                 % 扫描角度间隔
theta_scan = -90:dtheta:90;                 % 俯仰扫描角度,-90~90
phi_scan = -90:dphi:90;                     % 方位扫描角度，-90~90
theta_len = length(theta_scan);
phi_len = length(phi_scan);
beam = zeros(theta_len, phi_len);           % 初始化波束
for  r= 1:1:range
    for j = 1:1:phi_len
       theta = 45;
        phi = phi_scan(j);
        Fx = exp(1i*2*pi./lamda).*(X*sin(theta)*cos(phi)-(X.^2.*(1-(cos(phi)).^2.*(sin(theta)).^2))./2*r);%x axis
        Fy = exp(1i*2*pi./lamda).*(Y*cos(phi)-(Y.^2.*(sin(theta)).^2)./2*r);%y axis
        Fxy = kron(Fx,Fy); 
       beam(r,j) = abs(((axy.')'*(Fxy.')));

 %for  i = 1:1:theta_len
  %  for j = 1:1:phi_len
    %    theta = theta_scan(i);
     %   phi = phi_scan(j);
     %   Fx = exp(1i*2*pi./lamda).*(X*sind(theta)*cosd(phi)-(X.^2.*(1-(cosd(phi)).^2.*(sind(theta)).^2))./2*range);%x axis
     %   Fy = exp(1i*2*pi./lamda).*(Y*cosd(phi)-(Y.^2.*(sind(theta)).^2)./2*range);%y axis
      %  Fxy = kron(Fx,Fy); 
      %  beam(i,j) = abs(((axy.')'*(Fxy.')));
        
    end
end
beam_db = 20*log10(beam/max(max(beam)));
 
%CR under CC design
yita=1/pi;
delta = [];
p=0:20:120;
for i=-0.3:0.1:0.4 %y
    for j=-0.67:0.1:0.75 %z
        delta = (2/3) * (atan( (j*i) / ((sin(theta0) * cos(phi0))*sqrt((sin(theta0) * cos(phi0))^2+j^2+i^2)))) + ((sin(theta0) * cos(phi0)) * i *j )/ (3*((sin(theta0) * cos(phi0))^2+i^2)*sqrt((sin(theta0) * cos(phi0))^2+i^2+j^2));
        delta = [delta delta];
    end
end

CRcc_25 = log2(db2mag(p)) + log2((1/(4*pi)) * db2mag(sum(delta)));
CRcc_23 = log2(1 + db2mag(p) * (1/(4*pi) .* db2mag(sum(delta))));


A=lamda^2/(4*pi);

h_buffer = [];
for i=-7:1:7 %y 
    for j=-7:1:7
       % r=62.5000*sqrt((i*0.001-0.3536)^2 + (j*0.001-sqrt(2)/2)^2 + 3/8);
       % H= 121.5 * sqrt( A * ((62.5000^3.*(sqrt(6)/4)^3+ (62.5000*(sqrt(6)/4))^2 * (62.5000*(sqrt(2)/2) - j * (lamda/2))^2) / (4*pi*r^5)) ) * exp(-1i*2*pi*r/lamda) ;
      r=10.*sqrt((i*0.00625-sin(pi/4)*sin(pi/6))^2 + (j*0.00625-sqrt(2)/2)^2 + 3/8);
      H= sqrt( A * ((10^3.*(sqrt(6)/4)^3+ (10^2*(sqrt(6)/4))^2 * (10*(sqrt(2)/2) - j * (lamda/2))^2) / (4*pi*r^5)) ) * exp(-1i*2*pi*r/lamda) ;
      h_buffer = [h_buffer H]; 
    end
end
h_norm = norm(h_buffer)^2;
CRcc_22=log2(1+db2mag(p).* h_norm );


h_buffer1 = [];
for i=-7:1:7 %y
    for j=-7:1:7
        r=5*sqrt((i*0.0125-sin(pi/4)*sin(-pi/6))^2 + (j*0.0125-sqrt(2)/2)^2 + 3/8);
        H=sqrt( A * ((5^3.*(sqrt(6)/4)^3+ (5*(sqrt(6)/4))^2 * (5*(sqrt(2)/2) - j * (lamda/2))^2) / (4*pi*r^5)) ) * exp(-1i*2*pi*r/lamda) ;
        h_buffer1 = [h_buffer1 H]; 
    end
end
h_norm1 = norm(h_buffer1)^2; %h_s


%CR under SC design
delta1 = [];
p=0:20:120;
for i=-0.26:0.01:0.45 %y
    for j=-0.7:0.1:0.7 %z
        delta1 = (2/3) * (atan( (j*i) / ((sin(theta0) * cos(phi0))*sqrt((sin(theta0) * cos(phi0))^2+j^2+i^2)))) + ((sin(theta0) * cos(phi0)) * i *j )/ (3*((sin(theta0) * cos(phi0))^2+i^2)*sqrt((sin(theta0) * cos(phi0))^2+i^2+j^2));
        delta1 = [delta1 delta1];
    end
end
rou= ((conj(h_norm) * h_norm1)^2) / (h_norm * h_norm1);
CRsc_46=log2(1+(db2mag(p).*rou).*db2mag(sum(delta1)));



figure;
plot(p, CRcc_25, 'LineStyle',':','Color','r', 'LineWidth',1.5); hold on;

plot(p, CRcc_23, 'LineStyle','-', 'Color','b', 'LineWidth',1.5, 'Marker','o');

plot(p, CRcc_22, 'LineStyle','-', 'Color','g', 'LineWidth',1.5, 'Marker','*');
plot(p, CRsc_46, 'LineStyle','-', 'LineWidth',1.5, 'Marker','*');

legend('Approximation', 'Eq. (23)', 'Eq. (22)','Eq.(46)');
xlabel('p [dB]');
ylabel('CR bps/Hz');
grid on;



%SR under CC design

p=40:20:140;
L=4;

SRcc_41 = (1/L).*(log2(db2mag(p)) + 2*(log2( sqrt(L * 5.2143e-04)/(4*pi) * db2mag(sum(delta) + 21) )) );% does not follow?

SRcc_33 = (1/L).*log2( 1 + (db2mag(p).* (L).* rou ./ (16*pi^2)) * db2mag(sum(delta) + 21)^2);
SRcc_32 = (1/L).*log2( 1 + db2mag(p).* L.* db2mag(h_norm1)^2 .*rou);





%SR under SC design
SRsc_39 = (1/L).*log2( 1 + db2mag(p).* L.* db2mag(h_norm1)^2 .* 5.2143e-04);

SRcc_40 = (1/L).*log2( 1 + 5.2143e-04 .* (db2mag(p).* (L) ./ (16*pi*pi)) * db2mag(sum(delta1) + 21)^2 );


figure;
plot(p, SRcc_41, 'LineStyle',':','Color','k', 'LineWidth',1.5); hold on;

plot(p, SRcc_33, 'LineStyle','-', 'Color','r', 'LineWidth',1.5, 'Marker','o');

plot(p, SRsc_39, 'LineStyle','-', 'Color','b','LineWidth',1.5, 'Marker','*');

plot(p, SRcc_32, 'LineStyle','-',  'LineWidth',1.5, 'Marker','d');

plot(p, SRcc_40, 'LineStyle',':',  'LineWidth',1.5, 'Marker','d');


legend('Approximation', 'Eq. (33)', 'Eq. (39)', 'Eq. (32)', 'Eq. (40)');
xlabel('p [dB]');
ylabel('SR bps/Hz');
grid on;




%CR under SC design
rou=(conj(h_norm).'.*h_norm).^2./h_norm.^2;

CRsc_46=log2(1+(db2mag(p).*rou.*db2mag(yita/4*pi)).*db2mag(sum(delta)));



%CR-N 曲线
h2_buffer = [];
%N from 0 to 10^7, N=i*j, h relate to i,j; CR relate to h
for N=100:50:500
  for i1=(1-N)/2:10:(N-1)/2 %Ny
    for j1=(1-N)/2:10:(N-1)/2%Nz
        r=5*sqrt((i1*0.0125-0.3536)^2 + (j1*0.0125-sqrt(2)/2)^2 + 3/8);
        H1=sqrt( A * ((5^3.*(sqrt(6)/4)^3+ (5*(sqrt(6)/4))^2 * (5*(sqrt(2)/2) - j1 * (lamda/2))^2) / (4*pi*r^5)) ) * exp(-1i*2*pi*r/lamda) ;
        h2_buffer = [h2_buffer H1]; 
      
    end
  end
end
h2_norm_cal = [];
for i2=100:50:500
    h2_norm_cal = [h2_norm_cal norm(h2_buffer(i2:i2+6))];
end

CRcc_22_inf =log2(1 + db2mag(90).* h2_norm_cal );

figure;
plot(100:50:500,CRcc_22_inf,'LineStyle', '-', 'Color','r');

R_56_buffer = [];
R_57_buffer = [];
for tau = 0:0.1:1
    R_56 = (1/L) * log2(1 + db2mag(90) * L * h_norm1* ((tau^2 * rou * h_norm * h_norm1...
        + (1-tau)^2 * db2mag(h_norm1^2) + 2 * tau * (1-tau) * db2mag(h_norm) * db2mag(h_norm1^(1.5)))) / (tau^2 * db2mag(h_norm) ...
        + (1-tau)^2 * db2mag(h_norm1) + 2*tau*(1-tau)*db2mag(h_norm * h_norm1)));

    R_56_buffer = [R_56_buffer R_56];
    
    R_57 = log2(1 + db2mag(90)  * ((tau^2 * db2mag(h_norm^2) ...
        + (1-tau)^2 * rou * h_norm1*h_norm + 2 * tau * (1-tau) * db2mag(h_norm1) * db2mag(h_norm)^(1.5)))...
        / (tau^2 * db2mag(h_norm) + (1-tau)^2 * db2mag(h_norm1) + 2*tau*(1-tau)*db2mag(h_norm.^(0.5)) * db2mag(h_norm1.^(0.5))));

    R_57_buffer = [R_57_buffer R_57];
end
figure;
plot(R_57_buffer,R_56_buffer,'LineStyle', ':', 'Color', 'g','LineWidth',1.5, 'Marker','.');
legend('Downlink rate region');
xlabel('CR[bps/Hz]');
ylabel('SR bps/Hz');
line([0 0],[0 1.76]);
line([0 15], [0 0]);
hold on
grid on;

%R_58_buffer = [];
%R_59_buffer = [];
%for k=0.1:0.1:1
  %  for l=0:0.1:1
     %   R_58=k/L.*log2(1+db2mag((l/k).*p.*L).*h_norm1.^2.*rou);
     %   R_58_buffer = [R_58_buffer R_58];
     %   R_59=(1-k).*log2(1+db2mag(((1-l)./(1-k)).*p).*h_norm);
     %   R_59_buffer = [R_59_buffer R_59];

   % end
%end
%figure;
%plot(R_59_buffer,R_58_buffer,'LineStyle', ':', 'Color', 'g','LineWidth',1.5, 'Marker','.');

%grid on;

%CC design rc=2r
h3_buffer = [];
h4_buffer = [];
ha_buffer=[];
for ri=1:1:25
    for i=2.5:1:3 % Ny= 1001
      for j=-3:1:3 %Nz=1001
        r=2*ri*sqrt((i*(0.0625./2*ri)-0.3536)^2 + (j*(0.0625./2*ri)-sqrt(2)/2)^2 + 3/8);   
        h_USW=sqrt(A/(4*pi*(ri^2))) .*exp(-1i*(2*pi/lamda)*r);
        h_NUSW=sqrt(A/(4*pi*(r^2))).*exp(-1i*2*pi.*r/lamda); 
        h_acc=sqrt( A * (((2*ri)^3.*(sqrt(6)/4)^3+ ((2*ri)^2*(sqrt(6)/4))^2 * (2*ri*(sqrt(2)/2) - j * (lamda/2))^2) / (4*pi*(2*ri)^5)) ) * exp(-1i*2*pi*r/lamda);
        h4_buffer=[h4_buffer h_NUSW];
        h3_buffer = [h3_buffer h_USW]; 
        ha_buffer=[ha_buffer h_acc];

      end
    end
end

h3_norm_cal = [];
for ii=1:7:175
    h3_norm_cal = [h3_norm_cal norm(h3_buffer(ii:ii+6))];
end


h4_norm_cal = [];
for ii=1:7:175
    h4_norm_cal = [h4_norm_cal norm(h4_buffer(ii:ii+6))];
end

ha_norm_cal = [];
for ii=1:7:175
    ha_norm_cal = [ha_norm_cal norm(ha_buffer(ii:ii+6))];
end

CR22_NUSW=log2(1+db2mag(85).*(h4_norm_cal.^2));
CR22_USW=log2(1+db2mag(85).*(h3_norm_cal.^2));
CR22_acc=log2(1+db2mag(85).*(ha_norm_cal.^2));


figure;
plot(1:1:25, CR22_USW, 'LineStyle',':','LineWidth',1.5); hold on;
plot(1:1:25, CR22_NUSW, 'LineStyle','-','LineWidth',1.5); hold on;
plot(1:1:25, CR22_acc, 'LineStyle','-','LineWidth',1.5); hold off

legend('USW/UPW', 'NUSW', 'accurate');

%SC under r
%SC design rs=r
h5_buffer = [];
h6_buffer = [];
hacc_buffer=[];
for ri_s=2:2:50
    for i=2.5:1:3 % Ny= 1001
      for j=-3:1:3 %Nz=1001
        r=ri_s*sqrt((i*(0.0625./ri_s)-0.3536)^2 + (j*(0.0625./2*ri_s)-sqrt(2)/2)^2 + 3/8);   
        h_USW_s=sqrt(A/(4*pi*(ri_s^2))) .*exp(-1i*(2*pi/lamda)*r);
        h_NUSW_s=sqrt(A/(4*pi*(r^2))).*exp(-1i*2*pi.*r/lamda); 
        h_acc_s=sqrt( A * ((ri_s^3.*(sqrt(6)/4)^3+ (ri_s^2*(sqrt(6)/4))^2 * (ri_s*(sqrt(2)/2) - j * (lamda/2))^2) / (4*pi*ri_s^5)) ) * exp(-1i*2*pi*r/lamda);
        h6_buffer=[h6_buffer h_NUSW_s];
        h5_buffer = [h5_buffer h_USW_s]; 
        hacc_buffer=[hacc_buffer h_acc];

      end
    end
end
h5_norm_cal = [];
for ii=1:7:175
    h5_norm_cal = [h5_norm_cal norm(h5_buffer(ii:ii+6))];
end


h6_norm_cal = [];
for ii=1:7:175
    h6_norm_cal = [h6_norm_cal norm(h6_buffer(ii:ii+6))];
end

hacc_norm_cal = [];
for ii=1:7:175
    hacc_norm_cal = [hacc_norm_cal norm(hacc_buffer(ii:ii+6))];
end

CR22_NUSW_s=log2(1+db2mag(85).*(h6_norm_cal.^2));
CR22_USW_s=log2(1+db2mag(85).*(h5_norm_cal.^2));
CR22_acc_s=log2(1+db2mag(85).*(hacc_norm_cal.^2));

figure;
plot(2:2:50, CR22_USW_s, 'LineStyle',':','LineWidth',1.5); hold on;
plot(2:2:50, CR22_NUSW_s, 'LineStyle','-','LineWidth',1.5); hold on;
plot(2:2:50, CR22_acc_s, 'LineStyle','-','LineWidth',1.5); hold off

legend('USW/UPW', 'NUSW', 'accurate');
%    CRsc_46_N=log2(1+(db2mag(p1).*rou.*db2mag(yita/4*pi)).*db2mag(sum(delta2)));
%    figure;
%plot(p, CRsc_46_N, 'LineStyle',':','LineWidth',1.5); hold on;

%CCF vs N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%uplink CR under SC&CC design

pc=25:5:105;%ps=85dB
CRsc_69 = log2(1 + 4*db2mag(pc).* db2mag(sum(delta))./ ( 16*pi^2 + sum(delta1).^2.*db2mag(85)) ) ; % does not follow
CRsc_68 = log2(1 + db2mag(pc)*db2mag(h_norm) ./(1+db2mag(85).*db2mag(h_norm1).^2) );% too slow?
CRcc_up = log2(1 + db2mag(pc) * h_norm1);
CRsc_70 = log2(db2mag(pc))+log2( 4*pi.* sum(delta)./ ( 16*pi^2 +   sum(delta1).^2.*db2mag(85)) ) ; 



figure;
plot(pc, CRsc_68, 'LineStyle',':','Color','b', 'LineWidth',1.5, 'Marker','o'); hold on;
plot(pc, CRcc_up, 'LineStyle',':','Color','r', 'LineWidth',1.5,'Marker','*'); hold on;
plot(pc, CRsc_69, 'LineStyle',':','Color','k', 'LineWidth',1.5,'Marker','^'); hold on;
plot(pc, CRsc_70, 'LineStyle',':','Color','y', 'LineWidth',1.5); hold on;


legend('Eq. (68)', 'uplink-CRcc','CR-SC lower bound','high-SNR approximation' );
xlabel('pc [dB] with ps=85dB');
ylabel('CR bps/Hz');
grid on;


%uplink SR under CC design
ps=45:5:125;
SRcc_63 = (1/L).*log2(1+(db2mag(ps) .*L.*db2mag(h_norm1.^2).*(h_norm1-db2mag(60).*h_norm*h_norm1./(1+db2mag(60).*h_norm))) ) ;
SRcc_64 = (1/L).*(log2(1+db2mag(ps) *L.*yita^2.*sum(delta1).^2./(16*pi^2+4*pi.*yita*db2mag(60).*sum(delta) ) ) );
SRcc_65 = (1/L).*(log2(db2mag(ps)) +log2(L.*yita^2.*sum(delta).^2./(16*pi^2+4*pi.*yita*db2mag(60).*sum(delta) ) ) );%approx.



%uplink SR under SC design
SRsc_up=(1/L).*(log2(1+db2mag(ps) *L.*yita^2.*sum(delta1).^2./(16*pi^2 ) ) );



figure;
plot(ps, SRcc_63, 'LineStyle',':','Color','r', 'LineWidth',1.5,'Marker','^'); hold on;

plot(ps, SRcc_65, 'LineStyle','-', 'Color','k', 'LineWidth',1.5, 'Marker','o');
plot(ps, SRcc_64, 'LineStyle',':','Color','y', 'LineWidth',1.5); hold on;

plot(ps, SRsc_up, 'LineStyle','-', 'Color','b','LineWidth',1.5, 'Marker','*');

legend('ISAC-CC','high SNR_Approximation', 'ISAC CC-lower bound ', 'ISAC-SC');
xlabel('ps [dB] with pc=60dB');
ylabel('SR bps/Hz');
grid on;

 
ax = exp(1i*X*alpha_x);                 % channel vector in x-axis
ay = exp(1i*Y*alpha_y);                 % channel vector in x-axis
axy = kron(ax,ay);                      % Rectangular Formation Orientation Vector
 
 
dtheta = 0.2;
dphi = 0.2;                                 % intervals
theta_scan = -90:dtheta:90;                 % range of elevation,-90~90
phi_scan = -90:dphi:90;                     % range of azumith，-90~90
theta_len = length(theta_scan);
phi_len = length(phi_scan);
beam = zeros(theta_len, phi_len);           % utilizing beamforming
for i = 1:1:theta_len
    for j = 1:1:phi_len
        theta = theta_scan(i);
        phi = phi_scan(j);
        Fx = exp(1i*X*2*pi*sind(theta)*cosd(phi));
        Fy = exp(1i*Y*2*pi*sind(phi));
        Fxy = kron(Fx,Fy); 
        beam(i,j) = abs(((axy.')'*(Fxy.')));
    end
end
beam_db = 20*log10(beam/max(max(beam)));
 
figure;
mesh(phi_scan, theta_scan, beam_db);
title('Radiation pattern');
xlabel('Elevation');ylabel('azimuth');zlabel('amplitude(dB)');
axis([-100 100 -100 100 -80 10]);
 
figure;
imagesc(theta_scan,phi_scan,beam_db);
colorbar;axis tight;shading interp;
xlabel('Elevation');ylabel('azimuth');zlabel('amplitude(dB)');
title('rectangular surface array orientation diagram');

%beamforming
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
    lambda = 1;%lambda is unit 
    d = lambda / 2;
    theta = (0 : 0.001 : 2 * pi);
    V_theta = exp(-1i * (n - (N - 1) / 2)' * 2 * pi / lambda * cos(theta) * d);
    G_dB = 20 * log10(abs(ones(1, N) * V_theta / N));
end
 