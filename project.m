clc;clear;

% Load dataset
X=load('400101078.mat'); 

%struct to matrix
X=X.recieve;

% Parameters
K = size(X, 2); % Number of samples
N = size(X, 1); % Number of antennas
M = 2; % Number of radio sources

% ESPRIT algorithm

% finding correlation matrix
Rxx=X*(X.')./K; 

%eigen decompsition of Rxx
[Q, lambda] = eig(Rxx);
 
%sort to access the largest eighen values
[~, idx] = sort(diag(lambda), 'descend');

%partitioning using eighen values
Q_s = Q(:, idx(1:M)); 
Q_0=Q_s(1:end-1, :);
Q_1=Q_s(2:end, :);

%estimating psi using LS algorithm
psi = Q_0\Q_1;

%finding eighen values of psi
phi=eig(psi);

%finding estimated angles
theta = asind(angle(phi)/(pi));

% Display estimated angles
disp('Estimated angles of arrival:');
disp(theta);
