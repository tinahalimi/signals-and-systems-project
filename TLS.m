clc; clear;

% Load dataset
X = load('400101078.mat');
X = X.recieve;

% Parameters
K = size(X, 2); % Number of samples
N = size(X, 1); % Number of antennas
M = 2; % Number of radio sources

% ESPRIT algorithm with TLS
Rxx = X * (X.') / K; % Correlation matrix
[Q, lambda] = eig(Rxx); % Eigen decomposition of Rxx

[~, idx] = sort(diag(lambda), 'descend');

Q_s = Q(:, idx(1:M)); % Partitioning using eigenvalues
Q_0 = Q_s(1:N-1, :);
Q_1 = Q_s(2:N, :) - Q_s(1:N-1, :) * pinv(Q_s(1:N-1, :)) * Q_s(2:N, :); % Total least squares
psi = pinv(Q_0) * Q_1; % Estimating psi using TLS
phi = eig(psi); % Finding eigenvalues of psi
theta = asind(angle(phi)/(pi)); % Estimating angles

% Display estimated angles
disp('Estimated angles of arrival:');
disp(theta);