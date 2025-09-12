%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Traffic Matrix completion by Tensor Nuclear Norm
% Minus Forbenius Norm Minimization on Abilene Dataset 
%           30 Jan 2024, T.Miyata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;


close all;

addpath('mylib/'); %サブ関数 (HaLRTC由来)
addpath('utils/'); %テンソル関係の自作関数

%% Abilene読み込み-データセットはRoughanのウェブサイトから入手して前処理
load abilene_tm_2016.mat %2022年版
X=reshape(X,[144 2016]);



%% parameters
maxIter = 400;
%alpha = [1, 1, 1]; 0609　のとき
alpha = [1, 1, 1]; 
alpha = alpha / sum(alpha);
epsilon = 1e-30;
rho     = 1e-1;

RANK_SRSVD = 8; %----------SRSVDランク
ALPHA_GAIN = 1;

DIV=7; % データ点数が2016のとき

Prob_arr=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.92, 0.94, 0.96, 0.98];
%Prob_arr=[0.1,  0.5, 0.9, 0.98]; %時間がないとき

clear NMAE_H
clear NMAE_W
clear NMAE_S
clear NMAE_F

for i=1:numel(Prob_arr)
    % 乱数発生初期化
    rng('default');
    M_mask = rand(size(X)) >= Prob_arr(i);

    % source-destinationでテンソルにする
    T_org2 = mat2ten_SD(X);
    T_mask2 = logical(mat2ten_SD(M_mask));
    T2 = T_org2.*T_mask2;

    % 日ごとに分けてテンソルにする
    T_org = mat2ten_day(X, DIV);
    T_mask = logical(mat2ten_day(M_mask, DIV));
    T = T_org.*T_mask;

    M_org =X;
    M = M_org.*M_mask;

    %% HaLRTC
    [T_est_H, errList_H] = HaLRTC(T2, T_mask2, alpha,rho, maxIter,epsilon, T_org2);

    %% WLRTC
    [T_est_W, errList_W] = WLRTC_ADMM(T, T_mask, alpha*ALPHA_GAIN, rho,maxIter,epsilon,T_org); 

    %% TNNMFN-TC
    [T_est_F, errList_F] = TNNMFN_ADMM(T, T_mask, alpha*ALPHA_GAIN, rho,maxIter,epsilon,T_org);


    %% SRSVD
    [M_est_S err errobs]=ALS_part(M_mask, M, RANK_SRSVD, M_org, 0.1, 100);


    %% eval
    NMAE_H(i) = calc_NMAE(T_org, T_mask, T_est_H)
    NMAE_W(i) = calc_NMAE(T_org, T_mask, T_est_W)
    NMAE_S(i) = calc_NMAE(M_org, M_mask, M_est_S)
    NMAE_F(i) = calc_NMAE(T_org, T_mask, T_est_F)
end

%%
figure(1)
plot(Prob_arr, NMAE_S,'bs-','LineWidth',2,'MarkerSize',10)
hold on;
plot(Prob_arr, NMAE_H,'ko-','LineWidth',2,'MarkerSize',10)
plot(Prob_arr, NMAE_W,'v-','color', [0.3010 0.7450 0.9330],'LineWidth',2,'MarkerSize',10)
plot(Prob_arr, NMAE_F,'r*-','LineWidth',2,'MarkerSize',10)
xlim([0 1.0])
ylim([0 1.0])
legend('SRSVD', 'HaLRTC', 'WTTC-TS','TNFTC-TS(ours)','Location','northwest')
xlabel('Data missing probablity') % x-axis label
ylabel('NMAE') % y-axis label



pbaspect([1.3 1 1]);
ax = gca;
ax.FontName='Times New Roman';
ax.FontSize = 18;
saveas(gcf, sprintf('Abilene_NMAE_ANNET.eps'), 'epsc')
save result_Abilene_NMAE_ANNET.mat NMAE_S NMAE_H NMAE_W NMAE_F Prob_arr




