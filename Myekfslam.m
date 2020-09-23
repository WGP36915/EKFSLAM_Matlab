landmarks = [1 3 2 -1 -2; 1 2 3 4 1; 1 2 3 4 5];%每列代表一个landmark,最后一行是landmark的id
%landmarks = [1;1;1];
format compact;   %设置数据输出格式为压缩
screensize = get(0, 'ScreenSize')*0.6;  % ss = get(0,'ScreenSize')，返回行向量 ss = [ left, bottom, width, height ]
fig = figure('Position', [0 0 screensize(3), screensize(4)]);
plot(landmarks(1,:),landmarks(2,:),'b*') %画蓝色*代表Landmarks
hold on, axis equal, grid on
axis([-5 5 -1 5])%确定坐标系范围
    
    
%initialize states
global truePose ekfPose noisePose XX PX                          %全局变量
truePose = zeros(3,1);                                           %truePose是机器人的真实位姿
ekfPose = zeros(3,1);                                            %ekfPose是机器人带噪声位姿
noisePose = zeros(3,1);                                          %DR是机器人的带噪声位姿，不滤波
    
XX = zeros(3,1);                                                 %XX是EKF滤波中的状态变量   
PX = zeros(3);                                                   %PX代表协方差矩阵，初始化为3*3矩阵，即没有landmark时机器人自身位姿的协方差矩阵
    
linear_speed = 0.2;
angular_speed = 0.1;
dt = 0.1;
    
MAX_RANGE = 2;
lamda = 0.1;  %数据关联常数

%Control Noise
sigmaV = 0.2;
sigmaW = 1*pi/180.0;
Qt = [sigmaV^2 0 0; 0 sigmaV^2 0; 0 0 sigmaW^2];

%observation noise
sigmaR = 0.03;
sigmaB = (1.0*pi/180);
Rt = [sigmaR^2 0; 0 sigmaB^2];


associations = [];                                  %landmark的id集合

trueRobot = plot(0, 0, 'bo', 'erasemode', 'xor');
ekfRobot = plot(0, 0, 'ro', 'erasemode', 'xor');
noiseRobot = plot(0, 0, 'go', 'erasemode', 'xor');
observelandmark = plot(0, 0, 'r+', 'erasemode', 'xor');
legend('Landmark','Real Track','EKF-filtered Noise Track','Noise Track');
  while dt~=0
          %实际轨迹
          truePose(1) = truePose(1)+dt*linear_speed*cos( truePose(3) );
          truePose(2) = truePose(2)+dt*linear_speed*sin( truePose(3) );
          truePose(3) = pilimit(truePose(3)+angular_speed*dt);

          %带噪声的速度
          Vn = linear_speed+randn(1)*sigmaV;
          Wn = angular_speed+randn(1)*sigmaW;
          if(Wn~=0)
              nxUp = -(Vn/Wn)*sin(noisePose(3)) + (Vn/Wn)*sin(noisePose(3)+Wn*dt);
              nyUp = (Vn/Wn)*cos(noisePose(3)) - (Vn/Wn)*cos(noisePose(3)+Wn*dt);
          else
              nxUp = Vn*cos(noisePose(3));
              nyUp = Vn*sin(noisePose(3));  %这种运动模型只有在Wn等于
          end
          noisePose(1) = noisePose(1)+nxUp;
          noisePose(2) = noisePose(2)+nyUp;
          noisePose(3) = pilimit(noisePose(3)+Wn*dt);
         
          
          
         
          Gt = myJacobian(ekfPose(3), linear_speed, angular_speed, dt); 
          PX(1:3, 1:3) = Gt*PX(1:3, 1:3)*Gt'+Qt;                              %计算协方差矩阵
          if length(PX) > 3
              PX(1:3,4:end) = Gt*PX(1:3,4:end);
              PX(4:end,1:3) = PX(1:3,4:end)';
          end
          
          %预测
          if(Wn~=0)
              xUp = -(Vn/Wn)*sin(ekfPose(3)) + (Vn/Wn)*sin(ekfPose(3)+Wn*dt);
              yUp = (Vn/Wn)*cos(ekfPose(3)) - (Vn/Wn)*cos(ekfPose(3)+Wn*dt);
          else
              xUp = Vn*cos(ekfPose(3));
              yUp = Vn*sin(ekfPose(3));
          end
          
          ekfPose(1) = ekfPose(1)+xUp; %高斯噪声;
          ekfPose(2) = ekfPose(2)+yUp;
          ekfPose(3) = pilimit(ekfPose(3)+Wn*dt);
          
          XX(1) = ekfPose(1);
          XX(2) = ekfPose(2);
          XX(3) = ekfPose(3);
          
        
        %z = observed landmarks, zindex = it's id and relation to the association
        z = []; zindex = [];
        znew = true;
        for i=landmarks
            dist = sqrt((i(1) - truePose(1))^2 + (i(2) - truePose(2))^2);        %计算真实位置到路标的距离
            if dist <= MAX_RANGE                                           %距离小于MAX_RANGE表示看到了
                %All landmarks which can be observed described by range,bearing,id
                bearing = atan2( (i(2) - truePose(2)), (i(1) - truePose(1)) );

                l = [dist; pilimit(bearing - truePose(3)); i(3)];             %l存储路标相对机器人的位置，i(3)是landmark的id不变
                z = [z l];                                                 %加入机器人现在可以观测到的landmark集合中
                %Is this a new landmark? If not, get it's index in the map
                for y = 1:length(associations)
                    if i(3) == associations(y)
                        
                        zindex = [zindex y];
                        znew = false;
                    end
                end
                %Never seen this lm before, add it's ID to the association table
                if znew
                    associations = [associations i(3)];
                    %Increase the size of our mean and covar
                    muJ = [XX(1) + l(1)*cos(l(2) + XX(3)); XX(2) + l(1)*sin(l(2) + XX(3))]; %添加新的landmark在全局坐标系的位置作为新的状态
                    XX = [XX; muJ(1); muJ(2)];
                    covars = length(PX);
                    PX(covars + 1, covars + 1) = 1; 
                    PX(covars + 2, covars + 2) = 1;                        %协方差矩阵PX新增加的对角线元素设置为1(自己和自己的协方差)，其余默认位0(协方差为0，表示相互独立)
                end
            end
            znew = true;
        end
        
        
        %observation jacobian
        j = 0;
        for zi = 1:size(z,2)  %size(z,1)表示返回z的行数，size(z,2)表示返回列数
            %Calculated with current waypoint zloc
            zloc = z(:,zi);                                                %取矩阵z的第zi列，表示机器人现在观测到的第zi个landmark;
            for r = 1:length(associations) 
                if zloc(3) == associations(r)                              %zloc(3)是该landmark的id
                    j = r;
                end
            end
            if (j==0)
                disp('Cannot find landmark');
            end
             
            jind = j*2 + 2;                                                %2*j+2,2*j+3对应第j个landmark在状态XX中的 坐标值
            [observed, H] = observemodel(XX, jind); 
            
             set(observelandmark, 'xdata', XX(jind), 'ydata', XX(jind+1));  
             
            PHt = PX*H';
            S = H*PHt+Rt;                                                   %要特别注意是否可逆
            if(det(S)~=0)
               SInv = inv(S);
               K = PHt*SInv;            %卡尔曼增益
            else
               K = zeros(length(PX), 2); 
            end
            diff = [zloc(1) - observed(1); zloc(2) - observed(2)];                 %测量余差
            
            lam = diff'*S*diff; %数据关联值
            if(lam<lamda)
                XX = XX + K*[zloc(1) - observed(1); zloc(2) - observed(2)];            %更新状态
                PX = (eye(length(PX)) - K*H)*PX;
            end
            

            %有几个观测到的landmark就按顺序迭代几遍
            j=0;
        end
%              %画更新后的landmarks
%              idnum = (length(XX)-3)/2;
%              for id=1:idnum
%                   plot(XX(2*id+2), XX(2*id+3), 'r+', 'erasemode', 'xor');
%              end
           pause(dt/10);%延时
           
          ekfPose(1) = XX(1);
          ekfPose(2) = XX(2);
          ekfPose(2) = XX(2);
          
          
          set(trueRobot, 'xdata', truePose(1), 'ydata', truePose(2));
          plot(truePose(1), truePose(2), 'b');
          
          set(ekfRobot, 'xdata', ekfPose(1), 'ydata', ekfPose(2));
          plot(ekfPose(1), ekfPose(2), 'r');
          
          set(noiseRobot, 'xdata', noisePose(1), 'ydata', noisePose(2));
          plot(noisePose(1), noisePose(2), 'g');
          
%           legend('Landmark','Real Track','EKF-filtered Noise Track','Noise Track');
          pause(dt/10);%延时
  end
  
  