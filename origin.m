% 清空工作区和命令窗口
clear all;
clc;

% 参数设置
q = 0.1; % q的值

% 初始条件
c = @(x) 2*cos(pi * x); % w(x, 0)的初始条件
d = @(x) 2*sin(pi * x); % w_t(x, 0)的初始条件

% 时间参数
T = 30; % 模拟时间长度
dt = 0.01; % 时间步长，确保 dt = 0.01
t = 0:dt:T; % 时间网格
M = length(t); % 时间离散点数

% 空间参数
dx = 0.01; % 空间步长，确保 dx = 0.01
x = 0:dx:1; % 空间网格
N = length(x); % 空间离散点数

% 初始化变量
w = zeros(N, M); % 存储w(x, t)的矩阵
w_t = zeros(N, M); % 存储w_t(x, t)的矩阵

% 设置初始条件
w(:, 1) = c(x);
w_t(:, 1) = d(x);

% 计算 w(:, 2)
w(:, 2) = w(:, 1) + dt * w_t(:, 1);

% 模拟系统
for n = 2:M-1
    % 更新 w
    w(2:N-1, n+1) = -w(2:N-1, n-1) + w(3:N, n) + w(1:N-2, n);
    
    % 更新 w(1, n+1)
    w(1, n+1)=2*w(2, n)+ 2*dx*q*w(1, n)- w(1, n-1);
    
    % 更新 w(N, n+1)
    w(N, n+1)=2*w(N-1, n)-w(N, n-1);
    
end


% 绘制结果
figure;

mesh(t, x, w)
xlabel('t');
ylabel('x');
zlabel(['$w(x, t)$'],'Interpreter', 'latex');
title(['$w(x, t)$'],'Interpreter', 'latex');







