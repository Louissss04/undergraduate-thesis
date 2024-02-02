clear all;
clc;

% 参数设置
k1 = 15; 
k2 = 13; 

% 初始条件
c = @(x) cos(pi * x); 
d = @(x) sin(pi * x); 
z0 = k1 / k2; 

% 时间参数
T = 30; % 模拟时间长度
dt = 0.01; % 时间步长
t = 0:dt:T; % 时间网格
M = length(t); % 时间离散点数

% 空间参数
dx = 0.01; % 空间步长
x = 0:dx:1; % 空间网格
N = length(x); % 空间离散点数

% 初始化变量
w = zeros(N, M); % 存储w(x, t)的矩阵
wt = zeros(N, M); % 存储w_t(x, t)的矩阵
z = zeros(1, M); % 存储z(t)的矩阵

% 设置初始条件
w(:, 1) = c(x);
wt(:, 1) = d(x);
z(1) = z0;

% 更新 w(:, 2)
w(:, 2) = w(:, 1) + dt * wt(:, 1);
z(2) = z(1) + dt * (-k2 * z(1) + k2 * w(1, 1));

% 模拟系统
for n = 2:M-1
    % 更新w(x, t+dt)
    w(2:N-1, n+1) = -w(2:N-1, n-1) + w(3:N, n) + w(1:N-2, n);
    w(N, n+1) = w(N-1, n+1);
    w(1, n+1) = (w(2, n+1)+ dx * k2 *z(n))/(1+ dx* k1);
    
    % 更新w_t(x, t+dt)
    wt(:, n+1) = (w(:, n+1) - w(:, n)) / (dt);
    
    % 更新z(t+dt)
    z(n+1) = z(n) + dt * (-k2 * z(n) + k2 * w(1, n));
end

wt(:, 2) = (w(:, 3) - w(:, 1)) / (2*dt);%这一步是补上第二列wt

% 绘制结果
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


