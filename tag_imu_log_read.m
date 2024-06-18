clear all
close all
clc

g = 9.8155;

% load('imu_x_up.mat')
load('imu_x_up_kaf.mat')
a_meas_x_up = mean(tag_inf(:,2:4)/1e6*g);
w_meas_x_up = mean(tag_inf(:,5:7)/1e6*(pi/180));
m_meas_x_up = tag_inf(:,8:10)/1e6;
% load('imu_x_down.mat')
load('imu_x_down_kaf.mat')
a_meas_x_down = mean(tag_inf(:,2:4)/1e6*g);
w_meas_x_down = mean(tag_inf(:,5:7)/1e6*(pi/180));
m_meas_x_down = tag_inf(:,8:10)/1e6;
% load('imu_y_up.mat')
load('imu_y_up_kaf.mat')
a_meas_y_up = mean(tag_inf(:,2:4)/1e6*g);
w_meas_y_up = mean(tag_inf(:,5:7)/1e6*(pi/180));
m_meas_y_up = tag_inf(:,8:10)/1e6;
% load('imu_y_down.mat')
load('imu_y_down_kaf.mat')
a_meas_y_down = mean(tag_inf(:,2:4)/1e6*g);
w_meas_y_down = mean(tag_inf(:,5:7)/1e6*(pi/180));
m_meas_y_down = tag_inf(:,8:10)/1e6;
% load('imu_z_up.mat')
load('imu_z_up_kaf.mat')
a_meas_z_up = mean(tag_inf(:,2:4)/1e6*g);
w_meas_z_up = mean(tag_inf(:,5:7)/1e6*(pi/180));
m_meas_z_up = tag_inf(:,8:10)/1e6;
% load('imu_z_down.mat')
load('imu_z_down_kaf.mat')
a_meas_z_down = mean(tag_inf(:,2:4)/1e6*g);
w_meas_z_down = mean(tag_inf(:,5:7)/1e6*(pi/180));
m_meas_z_down = tag_inf(:,8:10)/1e6;

acc_bias_x = (a_meas_x_up(1) + a_meas_x_down(1))/2;
acc_bias_y = (a_meas_y_up(2) + a_meas_y_down(2))/2;
acc_bias_z = (a_meas_z_up(3) + a_meas_z_down(3))/2;
acc_bias = [acc_bias_x; acc_bias_y; acc_bias_z];

w_bias = mean([w_meas_x_up; w_meas_x_down; w_meas_y_up; w_meas_y_down; w_meas_z_up; w_meas_z_down])';


% load('imu_log_home_round.mat')

% ahrs_filter = IMU_AHRS_filter;
% ahrs_filter = ahrs_filter.config;

load('magcal_params_2.mat')
%load('magcal_params_1.mat')

clearvars -except magcal_params acc_bias w_bias

fid = fopen('imu_log_11042024165117','r');
k = 1;

while ~feof(fid)            
  currstr = fgetl(fid);
  if length(currstr) > 130
    tag_inf(k,:) = Decode_PDoA_tag(currstr);
    tag_str_log{k} = currstr;
    k = k + 1;
    clc
  end
end

fclose(fid);
g = 9.8155;
dw_unit = (1/499.2e6/128);
T_max = 2^40 * dw_unit;

for q = 1:length(tag_inf) 
    if q == 2500
       a = 1; 
    end
    a_meas(:,q) = tag_inf(q,2:4)/1e6*g;
    w_meas(:,q) = tag_inf(q,5:7)/1e6*(pi/180);
    m_meas(:,q) = tag_inf(q,8:10)/1e6;
    if q ~= 1
      Td(q,1) = (tag_inf(q,1) - tag_inf(q-1,1))*dw_unit;
      if Td(q,1) < 0
         Td(q,1) = tag_inf(q,1)*dw_unit + T_max - tag_inf(q-1,1)*dw_unit;  
      end
    else
      Td(q,1) = 0.02;  
    end
    mag_cal_x = m_meas(1,q)*1.0098 + 0.5085;
    mag_cal_y = m_meas(2,q)*0.9892 + 0.1433;
    mag_cal_z = m_meas(3,q)*1.0012 - 0.5797;
    m_meas_cal = [-mag_cal_x; mag_cal_y; mag_cal_z]; 
    axel(:,q)=a_meas(:,q) - acc_bias; 
    w(:,q)=w_meas(:,q)-w_bias;
    m(:,q)=m_meas_cal;
    dt(:,q)=Td(q,1);
end

% beta = -0.00142606;
beta = 0.025;
g_norm = [0;0;1];
h_norm = [0.3067189;0.06541019;0.9495499];
A = axel(:,1)/norm(axel(:,1));
M = m(:,1)/norm(m(:,1));
Mb = TRIAD(A,M);
Mr = TRIAD(g_norm,h_norm);
A_est = Mb*Mr';
Qt_1 = rotm2quat(A_est)';
Swt = zeros(1,4);
qet_dot = [0;0;0;0];
qt_stor = zeros(4,length(tag_inf));
qt_stor_2 = zeros(4,length(tag_inf));
Swct_stor = zeros(4,length(tag_inf));
Swbt_stor = zeros(4,length(tag_inf));
qt_norm_stor = zeros(1,length(tag_inf));
A_est_stor = zeros(3,3*length(tag_inf));
M_stor = zeros(3,3*length(tag_inf));
A_stor = zeros(3,length(tag_inf));
for n=1:length(tag_inf) 
    A(:,n) = axel(:,n)/norm(axel(:,n));
    M(:,n) = m(:,n)/norm(m(:,n));
    Swt(1,2)= w(1,n);
    Swt(1,3)= w(2,n);
    Swt(1,4)= w(3,n);
    dzeta = ((w(1,n)+w(2,n)+w(3,n))/3)*sqrt(3/4);
    qt2 = quat_sopr(Qt_1');
    Swet = 2*Quat_Mult(qt2,qet_dot');
    Swbt = dzeta*Swet*dt(1,n);
    Swct = Swt - Swbt';
    Swct(1,1) = 0;
    Swct_stor(:,n) = Swct;
    Swbt_stor(:,n) = Swbt';
    [qet_dot,qt]= madgwick(Qt_1',h_norm,A(:,n),M(:,n),Swct,beta,dt(1,n));
    qt_stor(1,n)=qt(1,1);
    qt_stor(2,n)=qt(2,1);
    qt_stor(3,n)=qt(3,1);
    qt_stor(4,n)=qt(4,1);
    qt_norm_stor(1,n) = norm(qt_stor(:,n));
    Qt_1(1,1)=qt(1,1);
    Qt_1(2,1)=qt(2,1);
    Qt_1(3,1)=qt(3,1);
    Qt_1(4,1)=qt(4,1);
end
for n=1:length(tag_inf) 
    A_1(:,n) = axel(:,n)/norm(axel(:,n));
    M_1(:,n) = m(:,n)/norm(m(:,n));
    Mb_2 = TRIAD(A_1(:,n),M_1(:,n));
    Mr_2 = TRIAD(g_norm,h_norm);
    M_stor(:,n) = Mb_2(:,1);
    M_stor(:,n+1) = Mb_2(:,1);
    M_stor(:,n+2) = Mb_2(:,1);
    A_est_2 = Mb_2*Mr_2';
    Quat_TRIAD = rotm2quat(A_est_2);
    qt_stor_2(1,n)=Quat_TRIAD(1,1);
    qt_stor_2(2,n)=Quat_TRIAD(1,2);
    qt_stor_2(3,n)=Quat_TRIAD(1,3);
    qt_stor_2(4,n)=Quat_TRIAD(1,4);
    A_est_stor(:,n) = A_est_2(:,1);
    A_est_stor(:,n+1) = A_est_2(:,2);
    A_est_stor(:,n+2) = A_est_2(:,3);
end
ahrs_filter = IMU_AHRS_filter;
ahrs_filter = ahrs_filter.config;
ahrs_filter.w_bias = w_bias;
for n=1:length(tag_inf)
    ahrs_filter = ahrs_filter.filter_update(a_meas(:,n)-acc_bias, w_meas(:,n), m_meas_cal ,Td(n,:));
    q_est_log(:,n) = ahrs_filter.q_est_out;
end
% Yaw_stor = rad2deg(Yaw_stor);
RPY = rad2deg(Quat2RPY(qt_stor));
figure
plot(RPY(1,:))
hold on
plot(RPY(2,:))
hold on
plot(RPY(3,:))
title('Roll Pitch Yaw Madgewick')
grid on

figure
plot(w_meas(1,:),'r')
hold on
plot(w_meas(2,:),'g')
hold on
plot(w_meas(3,:),'b')
hold on
plot(w_bias(1)*ones(1,length(a_meas)),'r')
hold on
plot(w_bias(2)*ones(1,length(a_meas)),'g')
hold on
plot(w_bias(3)*ones(1,length(a_meas)),'b')
title('Values and zero offsets')
grid on

figure
plot(qt_stor(1,:))
hold on
plot(qt_stor(2,:))
hold on
plot(qt_stor(3,:))
hold on
plot(qt_stor(4,:))
hold on
plot(qt_norm_stor(1,:))
title('Madgewick qaternion')
grid on

figure
plot(Swct_stor(1,:))
hold on
plot(Swct_stor(2,:))
hold on
plot(Swct_stor(3,:))
hold on
plot(Swct_stor(4,:))
hold on
title('Values of gyroscope including dzeta')
grid on

figure
plot(Swbt_stor(1,:))
hold on
plot(Swbt_stor(2,:))
hold on
plot(Swbt_stor(3,:))
hold on
plot(Swbt_stor(4,:))
hold on
title('Dzeta-driven displacements')
grid on

figure
plot(qt_stor_2(1,:))
hold on
plot(qt_stor_2(2,:))
hold on
plot(qt_stor_2(3,:))
hold on
plot(qt_stor_2(4,:))
title('TRIAD qaternion')
grid on

RPY_TRIAD = rad2deg(Quat2RPY(qt_stor_2));
figure
plot(RPY_TRIAD(1,:))
hold on
plot(RPY_TRIAD(2,:))
hold on
plot(RPY_TRIAD(3,:))
title('Roll Pitch Yaw TRIAD')
grid on

figure
plot(m(1,:))
hold on
plot(m(2,:))
hold on
plot(m(3,:))
grid on

RPY_val = rad2deg(Quat2RPY(q_est_log));
figure
plot(RPY_val(1,:))
hold on
plot(RPY_val(2,:))
hold on
plot(RPY_val(3,:))
title('Roll Pitch Yaw Valenti')
grid on