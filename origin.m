% ��չ������������
clear all;
clc;

% ��������
q = 0.1; % q��ֵ

% ��ʼ����
c = @(x) 2*cos(pi * x); % w(x, 0)�ĳ�ʼ����
d = @(x) 2*sin(pi * x); % w_t(x, 0)�ĳ�ʼ����

% ʱ�����
T = 30; % ģ��ʱ�䳤��
dt = 0.01; % ʱ�䲽����ȷ�� dt = 0.01
t = 0:dt:T; % ʱ������
M = length(t); % ʱ����ɢ����

% �ռ����
dx = 0.01; % �ռ䲽����ȷ�� dx = 0.01
x = 0:dx:1; % �ռ�����
N = length(x); % �ռ���ɢ����

% ��ʼ������
w = zeros(N, M); % �洢w(x, t)�ľ���
w_t = zeros(N, M); % �洢w_t(x, t)�ľ���

% ���ó�ʼ����
w(:, 1) = c(x);
w_t(:, 1) = d(x);

% ���� w(:, 2)
w(:, 2) = w(:, 1) + dt * w_t(:, 1);

% ģ��ϵͳ
for n = 2:M-1
    % ���� w
    w(2:N-1, n+1) = -w(2:N-1, n-1) + w(3:N, n) + w(1:N-2, n);
    
    % ���� w(1, n+1)
    w(1, n+1)=2*w(2, n)+ 2*dx*q*w(1, n)- w(1, n-1);
    
    % ���� w(N, n+1)
    w(N, n+1)=2*w(N-1, n)-w(N, n-1);
    
end


% ���ƽ��
figure;

mesh(t, x, w)
xlabel('t');
ylabel('x');
zlabel(['$w(x, t)$'],'Interpreter', 'latex');
title(['$w(x, t)$'],'Interpreter', 'latex');







