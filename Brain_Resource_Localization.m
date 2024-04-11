%% load data
clc;
clear;
close all;
load('ElecPosXYZ');
load("ElecPatch.mat");
load("Interictal.mat");

%% part a 
% Forward Matrix
ModelParams.R = [8 8.5 9.2];
ModelParams.Sigma = [3.3e-3 8.25e-5 3.3e-3];
ModelParams.Lambda = [.5979 .2037 .0237];
ModelParams.Mu = [.6342 .9364 1.0362];

Resolution = 1;
[LocMat, Gain] = ForwardModel_3shell(Resolution, ModelParams);

fig = figure('color',[1 1 1],'Renderer', 'painters', 'Position', [10 10 900 600]);
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'x'); 
axis('equal');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Dipoles positions');
saveas(fig,'Fig1.png');

%% part b,c
ElecNames = cell(1, 21);
ElecXYZ = zeros(3, 21);
for i = 1:21
    ElecNames{i} = ElecPos{1, i}.Name;
    ElecXYZ(:, i) = ElecPos{1, i}.XYZ';
end
ElecXYZ = ElecXYZ * 9.2;
fig = figure('color',[1 1 1],'Renderer', 'painters', 'Position', [10 10 900 600]);
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'x'); 
axis('equal');
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Dipoles/Electrodes positions');
hold on;
scatter3(ElecXYZ(1, :), ElecXYZ(2, :), ElecXYZ(3, :), 'o', 'filled');
for i = 1:21
    text(ElecXYZ(1, i)+0.5, ElecXYZ(2, i)+0.5, ElecXYZ(3, i)+0.5, ElecNames{i})
end
n = 1289;
e_q = LocMat(:, n) / norm(LocMat(:, n));  
line([LocMat(1,n),LocMat(1,n)+e_q(1)],[LocMat(2,n),LocMat(2,n)+e_q(2)],[LocMat(3,n),LocMat(3,n)+e_q(3)],'Color','green','linewidth',1.5);
scatter3(LocMat(1, n), LocMat(2, n), LocMat(3, n), 's', 'filled')
legend('Dipoles', 'Electrodes', 'Selected Dipole')
saveas(fig,'Fig2.png');

%% part d
G_select = Gain(:, 3*n-2:3*n);
e_q = LocMat(:, n) / norm(LocMat(:, n));
eq_superficial = e_q;
Q = e_q * Interictal(1, :);
M = G_select*Q;
disp_eeg(M, [], 256, ElecNames);

%% part e
max_amp = -inf;
max_idx = 0;
for i = 1:21
    if max(M(i, :) > max_amp)
        max_amp = max(M(i, :));
        max_idx = i;
    end
end

[~, locs] = findpeaks(M(max_idx, :), 'MinPeakDistance', 10, 'MinPeakHeight', 10);
mean_signal = zeros(21,1);
for i = 1:21
    for j = 1:length(locs)
        mean_signal(i) = mean_signal(i) + mean(M(i, locs(j)-3:locs(j)+3))/length(locs);
    end
end
fig = figure('color',[1 1 1],'Renderer', 'painters', 'Position', [10 10 900 600]);
scatter3(ElecXYZ(1, :), ElecXYZ(2, :), ElecXYZ(3, :), 36, mean_signal, 'o', 'filled');
colorbar;
colormap('jet');
for i = 1:21
    text(ElecXYZ(1, i)+0.5, ElecXYZ(2, i)+0.5, ElecXYZ(3, i)+0.5, ElecNames{i})
end
saveas(fig,'Fig4.png');

%% part f
Display_Potential_3D(9.2, mean_signal)
colormap('jet');

%% part g
alpha = 0.1;

% MNE
Q_MNE = Gain'*(Gain*Gain' + alpha*eye(21))^(-1)*mean_signal;

% WMNE
P = size(Gain, 2) / 3;
N = 21;
Omega = zeros(1, P);
for i = 1:P
    temp = 0;
    for j = 1:N
        g = Gain(j, 3*i-2:3*i);
        temp = temp + g*g';
    end
    Omega(i) = sqrt(temp);
end

W_WMNE = kron(diag(Omega), eye(3));
Q_WMNE = (W_WMNE'*W_WMNE)^(-1)*Gain'*(Gain*(W_WMNE'*W_WMNE)^(-1)*Gain' + alpha*eye(N))^(-1)*mean_signal;

% LORETA
d = 1;
A1 = zeros(P, P);
for alpha1 = 1:P
    for beta = 1:P
        A1(alpha1, beta) = (norm(LocMat(:, alpha1) - LocMat(:, beta)) == d)/6;
    end
end
A0 = diag(A1*ones(P, 1))^(-1)*A1;
A = kron(A0, eye(3));
B = 6/(d^2)*(A-eye(3*P));
W_LORETA = kron(diag(Omega), eye(3)) * (B' * B) * kron(diag(Omega), eye(3));
Q_LORETA = (W_LORETA'*W_LORETA)^(-1)*Gain'*(Gain*(W_LORETA'*W_LORETA)^(-1)*Gain' + alpha*eye(N))^(-1)*mean_signal;

% SLORETA
S_Q = Gain'*(Gain*Gain' + alpha*eye(N))^(-1)*Gain;
Q_SLORETA = zeros(3*P, 1);
for i = 1:P
    Q_SLORETA(3*i-2:3*i) = Q_MNE(3*i-2:3*i)'*S_Q(3*i-2:3*i, 3*i-2:3*i)^(-1)*Q_MNE(3*i-2:3*i);
end

%% part h
Q_MNE_max = -inf;
Q_MNE_max_idx = 0;
Q_WMNE_max = -inf;
Q_WMNE_max_idx = 0;
Q_LORETA_max = -inf;
Q_LORETA_max_idx = 0;
Q_SLORETA_max = -inf;
Q_SLORETA_max_idx = 0;
for i = 1:P
    if norm(Q_MNE(3*i-2:3*i)) > Q_MNE_max
        Q_MNE_max = norm(Q_MNE(3*i-2:3*i));
        Q_MNE_max_idx = i;
    end
    if norm(Q_WMNE(3*i-2:3*i)) > Q_WMNE_max
        Q_WMNE_max = norm(Q_WMNE(3*i-2:3*i));
        Q_WMNE_max_idx = i;
    end
    if norm(Q_LORETA(3*i-2:3*i)) > Q_LORETA_max
        Q_LORETA_max = norm(Q_LORETA(3*i-2:3*i));
        Q_LORETA_max_idx = i;
    end
    if norm(Q_SLORETA(3*i-2:3*i)) > Q_SLORETA_max
        Q_SLORETA_max = norm(Q_SLORETA(3*i-2:3*i));
        Q_SLORETA_max_idx = i;
    end
end

fig = figure('color',[1 1 1],'Renderer', 'painters', 'Position', [10 10 900 600]);
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'x'); 
axis('equal');
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Dipoles/Electrodes positions');
hold on;
scatter3(ElecXYZ(1, :), ElecXYZ(2, :), ElecXYZ(3, :), 'o', 'filled');
for i = 1:21
    text(ElecXYZ(1, i)+0.5, ElecXYZ(2, i)+0.5, ElecXYZ(3, i)+0.5, ElecNames{i})
end

scatter3(LocMat(1, n), LocMat(2, n), LocMat(3, n), 's', 'filled')
scatter3(LocMat(1, Q_MNE_max_idx), LocMat(2, Q_MNE_max_idx), LocMat(3, Q_MNE_max_idx), 's', 'filled')
scatter3(LocMat(1, Q_WMNE_max_idx), LocMat(2, Q_WMNE_max_idx), LocMat(3, Q_WMNE_max_idx), 's', 'filled')
scatter3(LocMat(1, Q_LORETA_max_idx), LocMat(2, Q_LORETA_max_idx), LocMat(3, Q_LORETA_max_idx), 's', 'filled')
scatter3(LocMat(1, Q_SLORETA_max_idx), LocMat(2, Q_SLORETA_max_idx), LocMat(3, Q_SLORETA_max_idx), 's', 'filled')
legend('Dipoles', 'Electrodes', 'Selected Dipole', 'Predict (MNE)', 'Predict (WMNE)', 'Predict (LORETA)', 'Predict (SLORETA)')
saveas(fig,'Fig6.png');

% location error
LocErr_MNE = mse(LocMat(:, n), LocMat(:, Q_MNE_max_idx)); %#ok<NASGU> 
LocErr_WMNE = mse(LocMat(:, n), LocMat(:, Q_WMNE_max_idx)); %#ok<NASGU> 
LocErr_LORETA = mse(LocMat(:, n), LocMat(:, Q_LORETA_max_idx)); %#ok<NASGU>
LocErr_SLORETA = mse(LocMat(:, n), LocMat(:, Q_SLORETA_max_idx)); %#ok<NASGU> 

% angle error
e_q_MNE_predict = Q_MNE(Q_MNE_max_idx-2:Q_MNE_max_idx) / norm(Q_MNE(Q_MNE_max_idx-2:Q_MNE_max_idx));
AngErr_MNE = acos((e_q'*e_q_MNE_predict) / (norm(e_q)*norm(e_q_MNE_predict)))*180/pi; %#ok<NASGU> 
e_q_WMNE_predict = Q_WMNE(Q_WMNE_max_idx-2:Q_WMNE_max_idx) / norm(Q_WMNE(Q_WMNE_max_idx-2:Q_WMNE_max_idx));
AngErr_WMNE = acos((e_q'*e_q_WMNE_predict) / (norm(e_q)*norm(e_q_WMNE_predict)))*180/pi; %#ok<NASGU> 
e_q_LORETA_predict = Q_LORETA(Q_LORETA_max_idx-2:Q_LORETA_max_idx) / norm(Q_LORETA(Q_LORETA_max_idx-2:Q_LORETA_max_idx));
AngErr_LORETA = acos((e_q'*e_q_LORETA_predict) / (norm(e_q)*norm(e_q_LORETA_predict)))*180/pi; %#ok<NASGU> 
e_q_SLORETA_predict = Q_SLORETA(Q_SLORETA_max_idx-2:Q_SLORETA_max_idx) / norm(Q_SLORETA(Q_SLORETA_max_idx-2:Q_SLORETA_max_idx));
AngErr_SLORETA = acos((e_q'*e_q_SLORETA_predict) / (norm(e_q)*norm(e_q_SLORETA_predict)))*180/pi;  %#ok<NASGU> 

%% part b
ElecNames = cell(1, 21);
ElecXYZ = zeros(3, 21);
for i = 1:21
    ElecNames{i} = ElecPos{1, i}.Name;
    ElecXYZ(:, i) = ElecPos{1, i}.XYZ';
end
ElecXYZ = ElecXYZ * 9.2;
fig = figure('color',[1 1 1],'Renderer', 'painters', 'Position', [10 10 900 600]);
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'x'); 
axis('equal');
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Dipoles/Electrodes positions');
hold on;
scatter3(ElecXYZ(1, :), ElecXYZ(2, :), ElecXYZ(3, :), 'o', 'filled');
for i = 1:21
    text(ElecXYZ(1, i)+0.5, ElecXYZ(2, i)+0.5, ElecXYZ(3, i)+0.5, ElecNames{i})
end
n = 853;
scatter3(LocMat(1, n), LocMat(2, n), LocMat(3, n), 's', 'filled') % dipole
legend('Dipoles', 'Electrodes', 'Selected Dipole')
saveas(fig,'Fig7.png');

%% part d
G_select = Gain(:, 3*n-2:3*n);
e_q = LocMat(:, n) / norm(LocMat(:, n));
e_q_deep = e_q;
Q = e_q * Interictal(1, :);
M = G_select*Q;
disp_eeg(M, [], 256, ElecNames);

%% part e
max_amp = -inf;
max_idx = 0;
for i = 1:21
    if max(M(i, :) > max_amp)
        max_amp = max(M(i, :));
        max_idx = i;
    end
end
[~, locs] = findpeaks(M(max_idx, :), 'MinPeakDistance', 10, 'MinPeakHeight', 10);
mean_signal = zeros(21,1);
for i = 1:21
    for j = 1:length(locs)
        mean_signal(i) = mean_signal(i) + mean(M(i, locs(j)-3:locs(j)+3))/length(locs);
    end
end

fig = figure('color',[1 1 1],'Renderer', 'painters', 'Position', [10 10 900 600]);
scatter3(ElecXYZ(1, :), ElecXYZ(2, :), ElecXYZ(3, :), 36, mean_signal, 'o', 'filled');
colorbar;
colormap('jet');
for i = 1:21
    text(ElecXYZ(1, i)+0.5, ElecXYZ(2, i)+0.5, ElecXYZ(3, i)+0.5, ElecNames{i})
end
hold on
scatter3(LocMat(1, n), LocMat(2, n), LocMat(3, n), 's', 'filled') % selected dipole
legend('Electrodes', 'Selected Dipole')
saveas(fig,'Fig9.png');

%% part f
Display_Potential_3D(9.2, mean_signal)
colormap('jet');

%% part g
alpha = 0.1;
% MNE
Q_MNE = Gain'*(Gain*Gain' + alpha*eye(21))^(-1)*mean_signal;

% WMNE
W_WMNE = kron(diag(Omega), eye(3));
Q_WMNE = (W_WMNE'*W_WMNE)^(-1)*Gain'*(Gain*(W_WMNE'*W_WMNE)^(-1)*Gain' + alpha*eye(N))^(-1)*mean_signal;

% LORETA
W_LORETA = kron(diag(Omega), eye(3)) * (B' * B) * kron(diag(Omega), eye(3));
Q_LORETA = (W_LORETA'*W_LORETA)^(-1)*Gain'*(Gain*(W_LORETA'*W_LORETA)^(-1)*Gain' + alpha*eye(N))^(-1)*mean_signal;

% SLORETA
Q_SLORETA = zeros(3*P, 1);

for i = 1:P
    Q_SLORETA(3*i-2:3*i) = Q_MNE(3*i-2:3*i)'*S_Q(3*i-2:3*i, 3*i-2:3*i)^(-1)*Q_MNE(3*i-2:3*i);
end

%% part h
Q_MNE_max = -inf;
Q_MNE_max_idx = 0;
Q_WMNE_max = -inf;
Q_WMNE_max_idx = 0;
Q_LORETA_max = -inf;
Q_LORETA_max_idx = 0;
Q_SLORETA_max = -inf;
Q_SLORETA_max_idx = 0;
for i = 1:P
    if norm(Q_MNE(3*i-2:3*i)) > Q_MNE_max
        Q_MNE_max = norm(Q_MNE(3*i-2:3*i));
        Q_MNE_max_idx = i;
    end
    if norm(Q_WMNE(3*i-2:3*i)) > Q_WMNE_max
        Q_WMNE_max = norm(Q_WMNE(3*i-2:3*i));
        Q_WMNE_max_idx = i;
    end
    if norm(Q_LORETA(3*i-2:3*i)) > Q_LORETA_max
        Q_LORETA_max = norm(Q_LORETA(3*i-2:3*i));
        Q_LORETA_max_idx = i;
    end
    if norm(Q_SLORETA(3*i-2:3*i)) > Q_SLORETA_max
        Q_SLORETA_max = norm(Q_SLORETA(3*i-2:3*i));
        Q_SLORETA_max_idx = i;
    end
end

fig = figure('color',[1 1 1],'Renderer', 'painters', 'Position', [10 10 900 600]);
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'x'); 
axis('equal');
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Dipoles/Electrodes Positions');
hold on;
scatter3(ElecXYZ(1, :), ElecXYZ(2, :), ElecXYZ(3, :), 'o', 'filled');
for i = 1:21
    text(ElecXYZ(1, i)+0.5, ElecXYZ(2, i)+0.5, ElecXYZ(3, i)+0.5, ElecNames{i})
end

scatter3(LocMat(1, n), LocMat(2, n), LocMat(3, n), 's', 'filled')
scatter3(LocMat(1, Q_MNE_max_idx), LocMat(2, Q_MNE_max_idx), LocMat(3, Q_MNE_max_idx), 's', 'filled')
scatter3(LocMat(1, Q_WMNE_max_idx), LocMat(2, Q_WMNE_max_idx), LocMat(3, Q_WMNE_max_idx), 's', 'filled')
scatter3(LocMat(1, Q_LORETA_max_idx), LocMat(2, Q_LORETA_max_idx), LocMat(3, Q_LORETA_max_idx), 's', 'filled')
scatter3(LocMat(1, Q_SLORETA_max_idx), LocMat(2, Q_SLORETA_max_idx), LocMat(3, Q_SLORETA_max_idx), 's', 'filled')
legend('Dipoles', 'Electrodes', 'Selected Dipole', 'Predict (MNE)', 'Predict (WMNE)', 'Predict (LORETA)', 'Predict (SLORETA)')
saveas(fig,'Fig11.png');

% location error
LocErr_MNE = mse(LocMat(:, n), LocMat(:, Q_MNE_max_idx));  
LocErr_WMNE = mse(LocMat(:, n), LocMat(:, Q_WMNE_max_idx));  
LocErr_LORETA = mse(LocMat(:, n), LocMat(:, Q_LORETA_max_idx)); 
LocErr_SLORETA = mse(LocMat(:, n), LocMat(:, Q_SLORETA_max_idx)); 

% angle error
e_q_MNE_predict = Q_MNE(Q_MNE_max_idx-2:Q_MNE_max_idx) / norm(Q_MNE(Q_MNE_max_idx-2:Q_MNE_max_idx));
AngErr_MNE = acos((e_q'*e_q_MNE_predict) / (norm(e_q)*norm(e_q_MNE_predict)))*180/pi;  
e_q_WMNE_predict = Q_WMNE(Q_WMNE_max_idx-2:Q_WMNE_max_idx) / norm(Q_WMNE(Q_WMNE_max_idx-2:Q_WMNE_max_idx));
AngErr_WMNE = acos((e_q'*e_q_WMNE_predict) / (norm(e_q)*norm(e_q_WMNE_predict)))*180/pi;  
e_q_LORETA_predict = Q_LORETA(Q_LORETA_max_idx-2:Q_LORETA_max_idx) / norm(Q_LORETA(Q_LORETA_max_idx-2:Q_LORETA_max_idx));
AngErr_LORETA = acos((e_q'*e_q_LORETA_predict) / (norm(e_q)*norm(e_q_LORETA_predict)))*180/pi;  
e_q_SLORETA_predict = Q_SLORETA(Q_SLORETA_max_idx-2:Q_SLORETA_max_idx) / norm(Q_SLORETA(Q_SLORETA_max_idx-2:Q_SLORETA_max_idx));
AngErr_SLORETA = acos((e_q'*e_q_SLORETA_predict) / (norm(e_q)*norm(e_q_SLORETA_predict)))*180/pi; 

%% part l Genetic Algorithm
x1 = [1289.0  -1.3882337679723866  -4.164700902593627  9.995282300223865];
q1 = x1(2:4);
e_q1 = q1 / norm(q1);
AngErr_GA = acos((eq_superficial'*e_q1') / (norm(e_q1')*norm(eq_superficial)))*180/pi; 
disp("AngErr_GA superficial dipole = "+AngErr_GA)

x2 = [853.0  6.23717949607848E-8  -5.373866997480836E-8  10.916849276991831];
q2 = x2(2:4);
e_q2 = q2/ norm(q2);
AngErr_GA_deep = acos((e_q_deep'*e_q2') / (norm(e_q2')*norm(e_q_deep)))*180/pi; 
disp("AngErr_GA deep dipole = "+AngErr_GA_deep)

%% part m
selected_dipoles = [1200:1204 1208:1212 1217:1221 1226:1230];

fig = figure('color',[1 1 1],'Renderer', 'painters', 'Position', [10 10 900 600]);
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'x'); 
axis('equal');
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Dipoles/Electrodes Positions');
hold on;
scatter3(ElecXYZ(1, :), ElecXYZ(2, :), ElecXYZ(3, :), 'o', 'filled');
for i = 1:21
    text(ElecXYZ(1, i)+0.5, ElecXYZ(2, i)+0.5, ElecXYZ(3, i)+0.5, ElecNames{i})
end
scatter3(LocMat(1, selected_dipoles), LocMat(2, selected_dipoles), LocMat(3, selected_dipoles),'MarkerFaceColor','green','MarkerEdgeColor','yellow')
e_q = zeros(3,length(selected_dipoles));

for i = 1:length(selected_dipoles)
    n = selected_dipoles(i);
    e_q(:,i) = LocMat(:,n)./norm(LocMat(:,n),2);
    n_eq = e_q(:,i);
    line([LocMat(1,n),LocMat(1,n)+n_eq(1)],[LocMat(2,n),LocMat(2,n)+n_eq(2)],[LocMat(3,n),LocMat(3,n)+n_eq(3)],'Color','yellow','linewidth',1.5);
end
legend('Dipoles', 'Electrodes', 'Selected dipole')
colormap('jet');
saveas(fig,'Fig12.png');

M = zeros(21,10240);
for i = 1:length(selected_dipoles)
    n = selected_dipoles(i);
    q = Interictal(i,:);
    Q = e_q(:,i)*q;
    M = M + Gain(:,3*n-2:3*n)*Q;
end

disp_eeg(M, [], 256, ElecNames);
title('Selected Dipoles');


max_amp = -inf;
max_idx = 0;
for i = 1:21
    if max(M(i, :) > max_amp)
        max_amp = max(M(i, :));
        max_idx = i;
    end
end

[~, locs] = findpeaks(M(max_idx, :), 'MinPeakDistance', 10, 'MinPeakHeight', 70);
mean_signal = zeros(21,1);
for i = 1:21
    for j = 1:length(locs)
        mean_signal(i) = mean_signal(i) + mean(M(i, locs(j)-3:locs(j)+3))/length(locs);
    end
end

fig = figure('color',[1 1 1],'Renderer', 'painters', 'Position', [10 10 900 600]);
scatter3(ElecXYZ(1, :), ElecXYZ(2, :), ElecXYZ(3, :), 36, mean_signal, 'o', 'filled');
colorbar;
for i = 1:21
    text(ElecXYZ(1, i)+0.5, ElecXYZ(2, i)+0.5, ElecXYZ(3, i)+0.5, ElecNames{i})
end
hold on
scatter3(LocMat(1, n), LocMat(2, n), LocMat(3, n), 's', 'filled')
legend('Electrodes', 'Selected Dipole');
colormap('jet');
saveas(fig,'Fig14.png');

Display_Potential_3D(9.2, mean_signal)
colormap('jet');

%% part n
alpha = 0.1;
% MNE
Q_MNE = Gain'*(Gain*Gain' + alpha*eye(21))^(-1)*mean_signal;

% WMNE
Q_WMNE = (W_WMNE'*W_WMNE)^(-1)*Gain'*(Gain*(W_WMNE'*W_WMNE)^(-1)*Gain' + alpha*eye(N))^(-1)*mean_signal;

% LORETA
W_LORETA = kron(diag(Omega), eye(3)) * (B' * B) * kron(diag(Omega), eye(3));
Q_LORETA = (W_LORETA'*W_LORETA)^(-1)*Gain'*(Gain*(W_LORETA'*W_LORETA)^(-1)*Gain' + alpha*eye(N))^(-1)*mean_signal;

% SLORETA
Q_SLORETA = zeros(3*P, 1);
for i = 1:P
    Q_SLORETA(3*i-2:3*i) = Q_MNE(3*i-2:3*i)'*S_Q(3*i-2:3*i, 3*i-2:3*i)^(-1)*Q_MNE(3*i-2:3*i);
end

Q_MNE1 = zeros(1317, 1);
Q_WMNE1 = zeros(1317, 1);
Q_LORETA1 = zeros(1317, 1);
Q_SLORETA1 = zeros(1317, 1);
for i = 1:P
    Q_MNE1(i) = norm(Q_MNE(3*i-2:3*i));
    Q_WMNE1(i) = norm(Q_WMNE(3*i-2:3*i));
    Q_LORETA1(i) = norm(Q_LORETA(3*i-2:3*i));
    Q_SLORETA1(i) = norm(Q_SLORETA(3*i-2:3*i));
end

%% part o
trueDipoles = zeros(1317, 1);
trueDipoles(selected_dipoles) = 1;
[X,Y] = perfcurve(trueDipoles,Q_MNE1,1);
fig = figure('color',[1 1 1],'Renderer', 'painters', 'Position', [10 10 900 600]);
plot(X, Y, 'Color', 'black'); 
xlabel('fpr'); 
ylabel('tpr');
hold on

[X,Y] = perfcurve(trueDipoles,Q_WMNE1,1); 
plot(X, Y, 'Color', 'blue');

[X,Y] = perfcurve(trueDipoles,Q_LORETA1,1); 
plot(X, Y, 'Color', 'red');

[X,Y] = perfcurve(trueDipoles,Q_SLORETA1,1); 
plot(X, Y, 'Color', 'green'); 
legend('MNE', 'WMNE', 'LORETA', 'sLORETA', 'interpreter', 'latex')
saveas(fig,'Fig16.png');

%% Functions
function t = disp_eeg(X,offset,feq,ElecName,titre)
    % function t = disp_eeg(X,offset,feq,ElecName,titre)
    %
    % inputs
    %     X: dynamics to display. (nbchannels x nbsamples) matrix
    %     offset: offset between channels (default max(abs(X)))
    %     feq: sapling frequency (default 1)
    %     ElecName: cell array of electrode labels (default {S1,S2,...})
    %     titre: title of the figure
    %
    % output
    %     t: time vector
    %
    % G. Birot 2010-02
    
    
    % Check arguments
    [N K] = size(X);
    
    if nargin < 4
        for n = 1:N
            ElecName{n}  = ['S',num2str(n)];
        end
        titre = [];
    end
    
    if nargin < 5
        titre = [];
    end
    
    if isempty(feq)
        feq = 1;
    end
    
    if isempty(ElecName)
        for n = 1:N
            ElecName{n}  = ['S',num2str(n)];
        end
    end
    
    if isempty(offset)
        offset = max(abs(X(:)));
    end
    
    
    % Build dynamic matrix with offset and time vector
    X = X + repmat(offset*(0:-1:-(N-1))',1,K);
    t = (1:K)/feq;
    graduations = offset*(0:-1:-(N-1))';
    shiftvec = N:-1:1;
    Ysup = max(X(1,:)) + offset;
    Yinf = min(X(end,:)) - offset;
    % YLabels = cell(N+2) ElecName(shiftvec)
    
    % Display
    figure1 = figure;
    % a1 = axes('YAxisLocation','right');
    a2 = axes('YTickLabel',ElecName(shiftvec),'YTick',graduations(shiftvec),'FontSize',7);
    ylim([Yinf Ysup]);
    box('on');
    grid('on')
    hold('all');
    plot(t,X');
    xlabel('Time (seconds)','FontSize',10);
    ylabel('Channels','FontSize',10);
    title(titre);
    hold off
end

function f = fitnessGA_fixed_eq(x)
    load('G.mat');
    load('LocMat.mat');
    load('mean_superficial.mat');
    n = x(1);
    q = x(2);
    e_q = LocMat(:,n) ./ norm(LocMat(:,n),2);
    GQ = G(:,3*n-2:3*n) * (q*e_q);
    
    f = norm(mean_signal - GQ,2);
end

function f = fitnessGA_notfixed_eq(x)
    load('G.mat');
    load('mean_deep.mat');
    n = x(1);
    q_eq = x(2:4);
    q_eq = reshape(q_eq,3,1);
    GQ = G(:,3*n-2:3*n) * q_eq;
    
    f = norm(mean_signal - GQ,2);
end