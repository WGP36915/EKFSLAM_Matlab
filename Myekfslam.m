landmarks = [1 3 2 -1 -2; 1 2 3 4 1; 1 2 3 4 5];%ÿ�д���һ��landmark,���һ����landmark��id
%landmarks = [1;1;1];
format compact;   %�������������ʽΪѹ��
screensize = get(0, 'ScreenSize')*0.6;  % ss = get(0,'ScreenSize')������������ ss = [ left, bottom, width, height ]
fig = figure('Position', [0 0 screensize(3), screensize(4)]);
plot(landmarks(1,:),landmarks(2,:),'b*') %����ɫ*����Landmarks
hold on, axis equal, grid on
axis([-5 5 -1 5])%ȷ������ϵ��Χ
    
    
%initialize states
global truePose ekfPose noisePose XX PX                          %ȫ�ֱ���
truePose = zeros(3,1);                                           %truePose�ǻ����˵���ʵλ��
ekfPose = zeros(3,1);                                            %ekfPose�ǻ����˴�����λ��
noisePose = zeros(3,1);                                          %DR�ǻ����˵Ĵ�����λ�ˣ����˲�
    
XX = zeros(3,1);                                                 %XX��EKF�˲��е�״̬����   
PX = zeros(3);                                                   %PX����Э������󣬳�ʼ��Ϊ3*3���󣬼�û��landmarkʱ����������λ�˵�Э�������
    
linear_speed = 0.2;
angular_speed = 0.1;
dt = 0.1;
    
MAX_RANGE = 2;
lamda = 0.1;  %���ݹ�������

%Control Noise
sigmaV = 0.2;
sigmaW = 1*pi/180.0;
Qt = [sigmaV^2 0 0; 0 sigmaV^2 0; 0 0 sigmaW^2];

%observation noise
sigmaR = 0.03;
sigmaB = (1.0*pi/180);
Rt = [sigmaR^2 0; 0 sigmaB^2];


associations = [];                                  %landmark��id����

trueRobot = plot(0, 0, 'bo', 'erasemode', 'xor');
ekfRobot = plot(0, 0, 'ro', 'erasemode', 'xor');
noiseRobot = plot(0, 0, 'go', 'erasemode', 'xor');
observelandmark = plot(0, 0, 'r+', 'erasemode', 'xor');
legend('Landmark','Real Track','EKF-filtered Noise Track','Noise Track');
  while dt~=0
          %ʵ�ʹ켣
          truePose(1) = truePose(1)+dt*linear_speed*cos( truePose(3) );
          truePose(2) = truePose(2)+dt*linear_speed*sin( truePose(3) );
          truePose(3) = pilimit(truePose(3)+angular_speed*dt);

          %���������ٶ�
          Vn = linear_speed+randn(1)*sigmaV;
          Wn = angular_speed+randn(1)*sigmaW;
          if(Wn~=0)
              nxUp = -(Vn/Wn)*sin(noisePose(3)) + (Vn/Wn)*sin(noisePose(3)+Wn*dt);
              nyUp = (Vn/Wn)*cos(noisePose(3)) - (Vn/Wn)*cos(noisePose(3)+Wn*dt);
          else
              nxUp = Vn*cos(noisePose(3));
              nyUp = Vn*sin(noisePose(3));  %�����˶�ģ��ֻ����Wn����
          end
          noisePose(1) = noisePose(1)+nxUp;
          noisePose(2) = noisePose(2)+nyUp;
          noisePose(3) = pilimit(noisePose(3)+Wn*dt);
         
          
          
         
          Gt = myJacobian(ekfPose(3), linear_speed, angular_speed, dt); 
          PX(1:3, 1:3) = Gt*PX(1:3, 1:3)*Gt'+Qt;                              %����Э�������
          if length(PX) > 3
              PX(1:3,4:end) = Gt*PX(1:3,4:end);
              PX(4:end,1:3) = PX(1:3,4:end)';
          end
          
          %Ԥ��
          if(Wn~=0)
              xUp = -(Vn/Wn)*sin(ekfPose(3)) + (Vn/Wn)*sin(ekfPose(3)+Wn*dt);
              yUp = (Vn/Wn)*cos(ekfPose(3)) - (Vn/Wn)*cos(ekfPose(3)+Wn*dt);
          else
              xUp = Vn*cos(ekfPose(3));
              yUp = Vn*sin(ekfPose(3));
          end
          
          ekfPose(1) = ekfPose(1)+xUp; %��˹����;
          ekfPose(2) = ekfPose(2)+yUp;
          ekfPose(3) = pilimit(ekfPose(3)+Wn*dt);
          
          XX(1) = ekfPose(1);
          XX(2) = ekfPose(2);
          XX(3) = ekfPose(3);
          
        
        %z = observed landmarks, zindex = it's id and relation to the association
        z = []; zindex = [];
        znew = true;
        for i=landmarks
            dist = sqrt((i(1) - truePose(1))^2 + (i(2) - truePose(2))^2);        %������ʵλ�õ�·��ľ���
            if dist <= MAX_RANGE                                           %����С��MAX_RANGE��ʾ������
                %All landmarks which can be observed described by range,bearing,id
                bearing = atan2( (i(2) - truePose(2)), (i(1) - truePose(1)) );

                l = [dist; pilimit(bearing - truePose(3)); i(3)];             %l�洢·����Ի����˵�λ�ã�i(3)��landmark��id����
                z = [z l];                                                 %������������ڿ��Թ۲⵽��landmark������
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
                    muJ = [XX(1) + l(1)*cos(l(2) + XX(3)); XX(2) + l(1)*sin(l(2) + XX(3))]; %����µ�landmark��ȫ������ϵ��λ����Ϊ�µ�״̬
                    XX = [XX; muJ(1); muJ(2)];
                    covars = length(PX);
                    PX(covars + 1, covars + 1) = 1; 
                    PX(covars + 2, covars + 2) = 1;                        %Э�������PX�����ӵĶԽ���Ԫ������Ϊ1(�Լ����Լ���Э����)������Ĭ��λ0(Э����Ϊ0����ʾ�໥����)
                end
            end
            znew = true;
        end
        
        
        %observation jacobian
        j = 0;
        for zi = 1:size(z,2)  %size(z,1)��ʾ����z��������size(z,2)��ʾ��������
            %Calculated with current waypoint zloc
            zloc = z(:,zi);                                                %ȡ����z�ĵ�zi�У���ʾ���������ڹ۲⵽�ĵ�zi��landmark;
            for r = 1:length(associations) 
                if zloc(3) == associations(r)                              %zloc(3)�Ǹ�landmark��id
                    j = r;
                end
            end
            if (j==0)
                disp('Cannot find landmark');
            end
             
            jind = j*2 + 2;                                                %2*j+2,2*j+3��Ӧ��j��landmark��״̬XX�е� ����ֵ
            [observed, H] = observemodel(XX, jind); 
            
             set(observelandmark, 'xdata', XX(jind), 'ydata', XX(jind+1));  
             
            PHt = PX*H';
            S = H*PHt+Rt;                                                   %Ҫ�ر�ע���Ƿ����
            if(det(S)~=0)
               SInv = inv(S);
               K = PHt*SInv;            %����������
            else
               K = zeros(length(PX), 2); 
            end
            diff = [zloc(1) - observed(1); zloc(2) - observed(2)];                 %�������
            
            lam = diff'*S*diff; %���ݹ���ֵ
            if(lam<lamda)
                XX = XX + K*[zloc(1) - observed(1); zloc(2) - observed(2)];            %����״̬
                PX = (eye(length(PX)) - K*H)*PX;
            end
            

            %�м����۲⵽��landmark�Ͱ�˳���������
            j=0;
        end
%              %�����º��landmarks
%              idnum = (length(XX)-3)/2;
%              for id=1:idnum
%                   plot(XX(2*id+2), XX(2*id+3), 'r+', 'erasemode', 'xor');
%              end
           pause(dt/10);%��ʱ
           
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
          pause(dt/10);%��ʱ
  end
  
  