clc;
clear;
%% ���������
g_NumOfPoints = 500;   % ����
g_ErrPointPart = 0.4;  % ����
g_NormDistrVar = 1;    % ��׼ƫ��
% ���������
theta = (rand(1) + 1) * pi/6;
R = ( rand([1 g_NumOfPoints]) - 0.5) * 100;
DIST = randn([1 g_NumOfPoints]) * g_NormDistrVar; %randn()������ֵΪ0������= 1����׼��= 1����̬�ֲ�������������ĺ���
Data = [cos(theta); sin(theta)] * R + [-sin(theta); cos(theta)] * DIST;

Data(:, 1:floor(g_ErrPointPart * g_NumOfPoints)) = 5*[10*rand(1), 10*rand(1); 10*rand(1), rand(1)] *(rand([2 floor(g_ErrPointPart * g_NumOfPoints)]) -randn([2 floor(g_ErrPointPart * g_NumOfPoints)])); %���������
% ... ��ʾMATLAB���ʽ��������һ�� 
% fix(x)����0ȡ��(Ҳ�������Ϊ���м�ȡ��), floor(x)������ȡ��, ceil(x)������ȡ��
plot(Data(1, :), Data(2, :), '.', 'Tag', 'DATA');
hold on;

%% RANSAC���
% ���ģ�ͳ�ʼ��
nSampLen = 2;                       %�趨ģ�������ݵĵ���
nIter = 50;                         %���ѭ������
dThreshold = 2;                     %�в���ֵ
nDataLen = size(Data, 2);           %���ݳ���
RANSAC_model = NaN;                 %����ȱʧģ��
RANSAC_mask = zeros([1 nDataLen]);  %ȫ0����1��ʾ����ģ�ͣ�0��ʾ������
nMaxInlyerCount = -1;
%% ��ѭ��
for i = 1:nIter 
    %  ������ѡȡ������ͬ�ĵ�
    SampleMask = zeros([1 nDataLen]);  
    while sum( SampleMask ) ~= nSampLen
        ind = ceil(nDataLen .* rand(1, nSampLen - sum(SampleMask))); %rand�����������ceil����������Ĵ�����Բ��
        SampleMask(ind) = 1;
    end    
    Sample = find( SampleMask );        %�ҳ�����Ԫ�ص�����ֵ��������ģ�͵ĵ�
    %% ����ģ��,�����ҷ���ģ�͵ĵ�
    ModelSet = feval(@TLS, Data(:, Sample));    %�������з���ģ�͵ĵ����С����
    for iModel = 1:size(ModelSet, 3) 
      CurModel = ModelSet(:, :, iModel);        %��ǰģ�Ͷ�Ӧ��ֱ�߲���   
      CurMask =( abs( CurModel * [Data; ones([1 size(Data, 2)])])< dThreshold);%��ֱ�߾���С����ֵ�ĵ����ģ��,���Ϊ1
      nCurInlyerCount = sum(CurMask);           %�������ֱ��ģ�͵ĵ�ĸ���
        %% ѡȡ���ģ��
        if nCurInlyerCount > nMaxInlyerCount    %����ģ�͵ĵ�������ģ�ͼ�Ϊ���ģ��
            nMaxInlyerCount = nCurInlyerCount;
            RANSAC_mask = CurMask;
            RANSAC_model = CurModel;
        end
    end
end
%% �����ģ�͵���Ͻ��
MinX=min(Data(1, :));
MaxX=max(Data(1, :));
MinX_Y=-(RANSAC_model(1).*MinX+RANSAC_model(3))./RANSAC_model(2);
MaxX_Y=-(RANSAC_model(1).*MaxX+RANSAC_model(3))./RANSAC_model(2);
plot([MinX MaxX],[MinX_Y MaxX_Y],'r');
title('RANSAC');
