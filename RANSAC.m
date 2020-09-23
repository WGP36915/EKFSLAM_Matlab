clc;
clear;
%% 随机点生成
g_NumOfPoints = 500;   % 点数
g_ErrPointPart = 0.4;  % 噪声
g_NormDistrVar = 1;    % 标准偏差
% 生成随机点
theta = (rand(1) + 1) * pi/6;
R = ( rand([1 g_NumOfPoints]) - 0.5) * 100;
DIST = randn([1 g_NumOfPoints]) * g_NormDistrVar; %randn()产生均值为0，方差= 1，标准差= 1的正态分布的随机数或矩阵的函数
Data = [cos(theta); sin(theta)] * R + [-sin(theta); cos(theta)] * DIST;

Data(:, 1:floor(g_ErrPointPart * g_NumOfPoints)) = 5*[10*rand(1), 10*rand(1); 10*rand(1), rand(1)] *(rand([2 floor(g_ErrPointPart * g_NumOfPoints)]) -randn([2 floor(g_ErrPointPart * g_NumOfPoints)])); %产生局外点
% ... 表示MATLAB表达式继续到下一行 
% fix(x)：向0取整(也可以理解为向中间取整), floor(x)：向左取整, ceil(x)：向右取整
plot(Data(1, :), Data(2, :), '.', 'Tag', 'DATA');
hold on;

%% RANSAC拟合
% 拟合模型初始化
nSampLen = 2;                       %设定模型所依据的点数
nIter = 50;                         %最大循环次数
dThreshold = 2;                     %残差阈值
nDataLen = size(Data, 2);           %数据长度
RANSAC_model = NaN;                 %跳过缺失模型
RANSAC_mask = zeros([1 nDataLen]);  %全0矩阵，1表示符合模型，0表示不符合
nMaxInlyerCount = -1;
%% 主循环
for i = 1:nIter 
    %  抽样，选取两个不同的点
    SampleMask = zeros([1 nDataLen]);  
    while sum( SampleMask ) ~= nSampLen
        ind = ceil(nDataLen .* rand(1, nSampLen - sum(SampleMask))); %rand产生随机数，ceil向离它最近的大整数圆整
        SampleMask(ind) = 1;
    end    
    Sample = find( SampleMask );        %找出非零元素的索引值，即建立模型的点
    %% 建立模型,并查找符合模型的点
    ModelSet = feval(@TLS, Data(:, Sample));    %计算所有符合模型的点的最小二乘
    for iModel = 1:size(ModelSet, 3) 
      CurModel = ModelSet(:, :, iModel);        %当前模型对应的直线参数   
      CurMask =( abs( CurModel * [Data; ones([1 size(Data, 2)])])< dThreshold);%到直线距离小于阈值的点符合模型,标记为1
      nCurInlyerCount = sum(CurMask);           %计算符合直线模型的点的个数
        %% 选取最佳模型
        if nCurInlyerCount > nMaxInlyerCount    %符合模型的点数最多的模型即为最佳模型
            nMaxInlyerCount = nCurInlyerCount;
            RANSAC_mask = CurMask;
            RANSAC_model = CurModel;
        end
    end
end
%% 画最佳模型的拟合结果
MinX=min(Data(1, :));
MaxX=max(Data(1, :));
MinX_Y=-(RANSAC_model(1).*MinX+RANSAC_model(3))./RANSAC_model(2);
MaxX_Y=-(RANSAC_model(1).*MaxX+RANSAC_model(3))./RANSAC_model(2);
plot([MinX MaxX],[MinX_Y MaxX_Y],'r');
title('RANSAC');
