% ��չ������������
clear all;
clc;

% ��������
k1 = 15; % k1��ֵ
k2 = 13; % k2��ֵ
c0 = 0.2; % c0��ֵ
c1 = 0.1; % c1��ֵ
q = 0.1; % q��ֵ

% ��ʼ����
c = @(x) 2*cos(pi * x); % w(x, 0)�ĳ�ʼ����
d = @(x) 2*sin(pi * x); % w_t(x, 0)�ĳ�ʼ����
chat = @(x) cos(pi * x); % w_hat(x, 0)�ĳ�ʼ����
dhat = @(x) sin(pi * x); % w_hat_t(x, 0)�ĳ�ʼ����
z0 = k1 / k2; % z(0)�ĳ�ʼ����

% ʱ�����
T = 300; % ģ��ʱ�䳤��
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
w_hat = zeros(N, M); % �洢w_hat(x, t)�ľ���
w_hat_t = zeros(N, M); % �洢w_hat_t(x, t)�ľ���
z = zeros(1, M); % �洢z(t)�ľ���

% ���ó�ʼ����
w(:, 1) = c(x);
w_t(:, 1) = d(x);
w_hat(:, 1) = chat(x);
w_hat_t(:, 1) = dhat(x);
z(1) = z0;

% ���� w(:, 2) �� wt(:, 2)
w(:, 2) = w(:, 1) + dt * w_t(:, 1);
w_hat(:, 2) = w_hat(:, 1) + dt * w_hat_t(:, 1);
z(2) = z(1) + dt * (-k2 * z(1) + k2 * w(1, 1) - k2 * w_hat(1, 1));

% ģ��ϵͳ
for n = 2:M-1
    % ���� w
    w(2:N-1, n+1) = -w(2:N-1, n-1) + w(3:N, n) + w(1:N-2, n);
    
    % ���� w_hat
    w_hat(2:N-1, n+1) = -w_hat(2:N-1, n-1) + w_hat(3:N, n) + w_hat(1:N-2, n);
    
    % ���� w(1, n+1)
    w(1, n+1)=2*w(2, n)+ 2*dx*q*w(1, n)- w(1, n-1);
    
    % ���� w_hat(1, n+1)
    w_hat(1, n+1)=(w_hat(2,n+1)+dx*(q+k1)*w(1,n+1)-dx*k2*z(n))/(1+dx*k1);
    
    
    % ���� w_hat(N, n+1)��w(N, n+1)
    S1=exp(q*(1-x)).*w_hat(:,n)';
    S2=exp(q*(1-x)).*w_hat_t(:,n)';
    w_hat(N, n+1)=(2*w_hat(N-1, n)-w_hat(N, n-1)-2*dx*(c0+q)*w_hat(N, n)-2*dx^2*(c0+q)*sum(S1,1:N-1)+c1*w_hat(N, n-1)-2*dx^2*c1*(c0+q)*sum(S2,1:N-1))/(1+c1);
    w(N, n+1)=w_hat(N, n+1);
    
    
    % ���� z
    z(n+1) = z(n) + dt * (-k2 * z(n) + k2 * (w(1, n) - w_hat(1, n)));
    
   
end

for n = 2:M-1
    % ���� w_t��w_hat_t
    w_t(:, n) = (w(:, n+1) - w(:, n-1)) / (2*dt);
    w_hat_t(:, n) = (w_hat(:, n+1) - w_hat(:, n-1)) / (2*dt);
end
    w_t(:, M) = (w(:, M) - w(:, M-1)) / (dt);
    w_hat_t(:, M) = (w_hat(:, M) - w_hat(:, M-1)) / (dt);

% ���ƽ��
figure;
subplot(2, 3, 1);
mesh(t, x, w)
xlabel('t');
ylabel('x');
zlabel(['$w(x, t)$'],'Interpreter', 'latex');
title(['$w(x, t)$'],'Interpreter', 'latex');

subplot(2, 3, 2);
mesh(t, x, w_hat);
xlabel('t');
ylabel('x');
zlabel(['$\hat{w}(x, t)$'],'Interpreter', 'latex');
title(['$\hat{w}(x, t)$'],'Interpreter', 'latex');

subplot(2, 3, 3);
mesh(t, x, w_t);
xlabel('t');
ylabel('x');
zlabel(['$w_t(x, t)$'],'Interpreter', 'latex');
title(['$w_t(x, t)$'],'Interpreter', 'latex');

subplot(2, 3, 4);
mesh(t, x, w_hat_t);
xlabel('t');
ylabel('x');
zlabel(['$\hat{w}_t(x, t)$'],'Interpreter', 'latex');
title(['$\hat{w}_t(x, t)$'],'Interpreter', 'latex');

subplot(2, 3, 5);
plot(t, z);
xlabel('t');
ylabel(['$z(t)$'],'Interpreter', 'latex');
title(['$z(t)$'],'Interpreter', 'latex');



