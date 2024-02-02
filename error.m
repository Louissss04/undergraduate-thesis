clear all;
clc;

% ��������
k1 = 15; 
k2 = 13; 

% ��ʼ����
c = @(x) cos(pi * x); 
d = @(x) sin(pi * x); 
z0 = k1 / k2; 

% ʱ�����
T = 30; % ģ��ʱ�䳤��
dt = 0.01; % ʱ�䲽��
t = 0:dt:T; % ʱ������
M = length(t); % ʱ����ɢ����

% �ռ����
dx = 0.01; % �ռ䲽��
x = 0:dx:1; % �ռ�����
N = length(x); % �ռ���ɢ����

% ��ʼ������
w = zeros(N, M); % �洢w(x, t)�ľ���
wt = zeros(N, M); % �洢w_t(x, t)�ľ���
z = zeros(1, M); % �洢z(t)�ľ���

% ���ó�ʼ����
w(:, 1) = c(x);
wt(:, 1) = d(x);
z(1) = z0;

% ���� w(:, 2)
w(:, 2) = w(:, 1) + dt * wt(:, 1);
z(2) = z(1) + dt * (-k2 * z(1) + k2 * w(1, 1));

% ģ��ϵͳ
for n = 2:M-1
    % ����w(x, t+dt)
    w(2:N-1, n+1) = -w(2:N-1, n-1) + w(3:N, n) + w(1:N-2, n);
    w(N, n+1) = w(N-1, n+1);
    w(1, n+1) = (w(2, n+1)+ dx * k2 *z(n))/(1+ dx* k1);
    
    % ����w_t(x, t+dt)
    wt(:, n+1) = (w(:, n+1) - w(:, n)) / (dt);
    
    % ����z(t+dt)
    z(n+1) = z(n) + dt * (-k2 * z(n) + k2 * w(1, n));
end

wt(:, 2) = (w(:, 3) - w(:, 1)) / (2*dt);%��һ���ǲ��ϵڶ���wt

% ���ƽ��
figure;

subplot(2, 2, 1);
mesh(t, x, w)
xlabel('t');
ylabel('x');
zlabel(['$w(x, t)$'],'Interpreter', 'latex');
title(['$w(x, t)$'],'Interpreter', 'latex');

subplot(2, 2, 2);
mesh(t, x, wt);
xlabel('t');
ylabel('x');
zlabel(['$w_t(x, t)$'],'Interpreter', 'latex');
title(['$w_t(x, t)$'],'Interpreter', 'latex');

subplot(2, 2, 3);
plot(t, z);
xlabel('t');
ylabel(['$z(t)$'],'Interpreter', 'latex');
title(['$z(t)$'],'Interpreter', 'latex');

subplot(2, 2, 4);
plot(t, wt(N, :));
xlabel('t');
ylabel(['$w_t(1, t)$'],'Interpreter', 'latex');
title(['$w_t(1, t)$'],'Interpreter', 'latex');


