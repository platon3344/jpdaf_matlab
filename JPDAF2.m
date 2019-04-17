%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 程序功能  :采用联合概率数据关联算法实现两个匀速运动目标的点迹与航迹的关联
% 输入变量  ：
%           target_position:  目标的初始位置
%           Pd:  检测概率
%           Pg： 正确量测落入跟踪门内的概率 
%           g_sigma： 门限
%           lambda： 杂波密度
%           target_delta=100： 观测噪声标准差
%           process_delta=2：  过程噪声标准差
%           T_target： 真实目标数（为2时有效）
%           - n：  采样次数 
%           - T： 采样间隔
%           -MC_number：仿真次数
% 输出变量  ：
%           无           
% 作者     ：zx
% 日期     ： 05-20-2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function JPDAF(target_position,Pd,Pg,g_sigma,lambda,target_delta,process_delta,T_target,n,T,MC_number)
tic

A = [1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];                   %状态转移矩阵
H = [1 0 0 0;0 0 1 0];                                   %观测矩阵  
R=[target_delta^2  0;0  target_delta^2];                 %观测协方差矩阵
Q=[process_delta^2 0;0 process_delta^2];                    %过程噪声协方差
G=[T^2/2 0;T 0;0 T^2/2;0 T];                             %过程噪声矩阵
P(:,:,1)=[target_delta^2 target_delta^2 0 0;target_delta^2 2*target_delta^2 0 0;0 0 target_delta^2 target_delta^2;0 0 target_delta^2 2*target_delta^2];
P(:,:,2)=[target_delta^2 target_delta^2 0 0;target_delta^2 2*target_delta^2 0 0;0 0 target_delta^2 target_delta^2;0 0 target_delta^2 2*target_delta^2];
                                                         %初始协方差矩阵

x_filter(:,:,1)=zeros(4,n);      %存储第一个目标的各时刻的滤波值
x_filter(:,:,2)=zeros(4,n);      %存储第二个目标的各时刻的滤波值
x_filter1(:,:,1)=zeros(4,n);     %MC_number次montle carlo仿真所得全部结果累加存储矩阵
x_filter1(:,:,2)=zeros(4,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%产生目标的实际位置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
data_measurement1=zeros(4,n,2);              	%data_measurement1实际位置矩阵
data_measurement1(:,1,1)=target_position(1,:)';
data_measurement1(:,1,2)=target_position(2,:)';
for tt=1:T_target
    for i=2:n
        data_measurement1(:,i,tt)=A*data_measurement1(:,i-1,tt)+G*sqrt(Q)*(randn(2,1));   %实际位置 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%主程序
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%产生路径
data_measurement=zeros(2,n,T_target);                  %观测存储矩阵,n采样次数
for tt=1:T_target
    for i=1:n
       data_measurement(:,i,tt)=H*data_measurement1(:,i,tt)+sqrt(R)*randn(2,1); 
%data_measurement观测值矩阵,传感器观测的位置	
    end
end
for M=1:MC_number
NOISE=[];
% 滤波开始
for t=1:n
    y1=[];
    for tt=1:T_target
        Noise=[];
        if t~=1
            x_predic(:,tt) = A*x_filter(:,t-1,tt); %用前一时刻的滤波值来预测当前的值 
        else
            x_predic(:,tt) = target_position(tt,:)'; %第一次采样用真实位置当预测值 
        end
        P_predic(:,:,tt) = A*P(:,:,tt)*A'+G*Q*G';
        Z_predic(:,tt) = H*x_predic(:,tt);
        S(:,:,tt) = H*P_predic(:,:,tt)*H'+ R;
        K(:,:,tt) = P_predic(:,:,tt)*H'*inv(S(:,:,tt));                         %增益
        ellipse_Volume(tt)=pi*g_sigma*sqrt(det(S(:,:,tt))); %计算椭球体积，这里算的是面积   
        number_returns(tt)=floor(10*ellipse_Volume(tt)*lambda+1); %错误回波数
        side(tt)=sqrt(10*ellipse_Volume(tt))/2; %求出正方行边长的二分之一
        Noise(1,:,tt)=data_measurement(1,t,tt)+side(tt)-2*rand(1,number_returns(tt))*side(tt); 
													%在预测值周围产生多余回波
        Noise(2,:,tt)=data_measurement(2,t,tt)+side(tt)-2*rand(1,number_returns(tt))*side(tt);
        NOISE=[NOISE Noise(:,:,tt)];
        %
        b=zeros(2,1);
        b=data_measurement(:,t,tt);
        y1=[y1 Noise(:,:,tt) b];            %将接收到的所有的回波存在y1中
    end
    y=[];
    d=[];
    Ohm=[];                                            %确认矩阵
    m=0;
    for j=1:number_returns(1)+number_returns(2)+2
        Ohm1=zeros(1,T_target+1);
        for tt=1:T_target
            d(:,tt)=y1(:,j)-Z_predic(:,tt);
            D(tt)=d(:,tt)'*inv(S(:,:,tt))*d(:,tt);
            if D(tt)<=g_sigma
                Ohm1(1,tt+1)=1;
            end
        end
        if sum(Ohm1(1,:))>0
           Ohm1(1,1)=1;
           Ohm=[Ohm;Ohm1];
           m=m+1;                        %记录总的有效回波个数
           y=[y y1(:,j)];              %把落入任意跟踪门中的所有回波放入y中
        end
    end
    
    %m=0表示无有效回波
         
    if m==0                                %无回波的情况
       x_filter(:,t,1)= x_predic(:,1);
       x_filter(:,t,2)=x_predic(:,2);
       P(:,:,1)=P_predic(:,:,1);                                           
       P(:,:,2)=P_predic(:,:,2);
    else
       if m*2-sum(Ohm(:))==0              %不同目标的有效回波没有关联的情况
           for tt=1:T_target
               yy=[];
               m=0;
               J=find(Ohm(:,tt+1)==1);
               m=length(J);
               if m==0
                   x_filter(:,t,tt)= x_predic(:,tt);
                   P(:,:,tt)=P_predic(:,:,tt); 
               else
                   for j=1:m
                       yy=[yy y(:,J(j))];
                   end
                    Bk=lambda*sqrt(2*pi*det(S(:,:,tt)))*(1-Pd*Pg)/Pd;      %算b0 
                    E=zeros(1,m);
                    belta=zeros(1,m);
                    for i=1:m
                        a=(yy(:,i)-Z_predic(:,tt))'*inv(S(:,:,tt))*(yy(:,i)-Z_predic(:,tt));
                        E(i)=E(i)+exp(-a/2);
                    end 
                    belta0=Bk/(Bk+sum(E));             %无回波时的关联概率
                    v=zeros(2,1);
                    v1=zeros(2,2);
                    for i=1:m
                        belta(i)=E(i)/(Bk+sum(E));                 %算关联概率
                        v=v+belta(i)*(yy(:,i)-Z_predic(:,tt));
                        v1=v1+belta(i)*(yy(:,i)-Z_predic(:,tt))*(yy(:,i)-Z_predic(:,tt))';
                    end
                    x_filter(:,t,tt)= x_predic(:,tt) + K(:,:,tt)*v;
                    %算协方差
                    Pc=(eye(4)-K(:,:,tt)*H)*P_predic(:,:,tt);
                    PP=K(:,:,tt)*(v1-v*v')*K(:,:,tt)';
                    P(:,:,tt)=belta0*P_predic(:,:,tt)+(1-belta0)*Pc+PP;
               end
           end
       else
                                             		%不同目标有回波互联的情况
           [Ohm_predict,Nry]=AResolution(Ohm);        %将确认矩阵拆分成互联矩阵
           c=0;                                     %归一化系数
           beta=zeros(m,T_target);                      %关联概率
           P_theta=ones(1,Nry);
           for i=1:Nry
               phi=sum(Ohm_predict(:,1,i));  
               for tt=1:T_target
                   delta(tt)=length(find(Ohm_predict(:,tt+1,i)==1));
                   P_theta(i)=P_theta(i)*(Pg*Pd)^(delta(tt))*(1-Pg*Pd)^(1-delta(tt));
               end
               [tau_x tau_y]=find(Ohm_predict(:,2:3,i)==1);
               for j=1:length(tau_x)
                   b=y(:,tau_x(j))-Z_predic(:,tau_y(j));
                   P_theta(i)=P_theta(i)*Pg^(-1)*(det(2*pi*S(:,:,tau_y(j))))^(-1/2)*exp(-(b'*inv(S(:,:,tau_y(j)))*b)/2);                   
               end
               P_theta(i)=P_theta(i)*lambda^phi;
               for j=1:length(tau_x)
                   beta(tau_x(j),tau_y(j))=beta(tau_x(j),tau_y(j))+P_theta(i);                   
               end          
               c=P_theta(i)+c;
           end
           if c==0
               for tt=1:T_target
                    Bk=lambda*sqrt(2*pi*det(S(:,:,tt)))*(1-Pd*Pg)/Pd;
                    E=(y-Z_predic(:,tt))'*inv(S(:,:,tt))*(y-Z_predic(:,tt));
                    beta(1:tt)=E/(Bk+sum(E));
               end
           else
               beta=beta/c;
           end

           for tt=1:T_target
                beta0=1-sum(beta(:,tt));
                v=zeros(2,1);
                v1=zeros(2,2);
                for i=1:m
                    v=v+beta(i,tt)*(y(:,i)-Z_predic(:,tt));
                    v1=v1+beta(i,tt)*(y(:,i)-Z_predic(:,tt))*(y(:,i)-Z_predic(:,tt))';
                end
                x_filter(:,t,tt)= x_predic(:,tt) + K(:,:,tt)*v;
                %算协方差
                Pc=(eye(4)-K(:,:,tt)*H)*P_predic(:,:,tt);
                PP=K(:,:,tt)*(v1-v*v')*K(:,:,tt)';
                P(:,:,tt)=beta0*P_predic(:,:,tt)+(1-beta0)*Pc+PP;
           end
       end
    end
end
x_filter1(:,:,1)=x_filter1(:,:,1)+x_filter(:,:,1);       %每次蒙特卡罗实验得到的滤波值作累加
x_filter1(:,:,2)=x_filter1(:,:,2)+x_filter(:,:,2);
end
toc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %画图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_filter(:,:,1)=x_filter1(:,:,1)/MC_number;                         %滤波值作平均
x_filter(:,:,2)=x_filter1(:,:,2)/MC_number;
figure                               %画出目标的估计位置，测量位置，杂波位置
color=['r-' 'r:' 'k-' 'k:'];
for tt=1:T_target
plot(x_filter(1,:,tt),x_filter(3,:,tt),color(2*tt-1:2*tt)),hold on
%plot(data_measurement1(1,:),data_measurement1(3,:),'-')
plot(data_measurement1(1,:,tt),data_measurement1(3,:,tt),color(2*tt+3:2*tt+4))
%axis([0 30 1 25])
end
plot(NOISE(1,:),NOISE(2,:),'.')                               %杂波位置
xlabel('x(m)'),ylabel('y(m)');
legend('目标1的估计位置','目标1的真实位置','目标2的估计位置','目标2的真实位置','杂波的位置',4)

for tt=1:T_target
    figure                                %画出目标X轴的估计位置，测量位置
    plot(1:n,x_filter(1,:,tt),'r-'),hold on
    plot(data_measurement1(1,:,tt),'k-')
    title(sprintf('目标%d',tt));
    xlabel('t(s)'),ylabel('x(m)');
    legend('X方向估计位置','X方向真实位置',4)
end

for tt=1:T_target
    figure                                %画出目标Y轴的估计位置，测量位置
    plot(1:n,x_filter(3,:,tt),'r-'),hold on
    plot(data_measurement1(3,:,tt),'k-')
    title(sprintf('目标%d',tt));
    xlabel('t(s)'),ylabel('y(m)');
    legend('Y方向估计位置','Y方向真实位置',4)
end

figure											%画出目标位置估计均方误差
color=['r:' 'b-'];
for tt=1:T_target
    a=zeros(1,n);                                   %位置均方误差 
    for j=1:n
a(1,j)=sqrt((x_filter(1,j,tt)-data_measurement1(1,j,tt))^2+(x_filter(3,j,tt)-data_measurement1(3,j,tt))^2);
    end
    plot(1:n,a(1,:),color(2*tt-1:2*tt));hold on;
    xlabel('t(s)'),ylabel('预测误差(m)');
end
legend('目标1的误差','目标2的误差',4)

figure											%画出目标位置估计均方误差
color=['r:' 'b-'];
for tt=1:T_target  
    a=zeros(1,n);                                    %速度均方误差
    for j=1:n
a(1,j)=sqrt((x_filter(2,j,tt)-data_measurement1(2,j,tt))^2+(x_filter(4,j,tt)-data_measurement1(4,j,tt))^2);
    end
    plot(1:n,a(1,:),color(2*tt-1:2*tt));hold on;
    xlabel('t(s)'),ylabel('预测误差(m/s)');
end
    legend('目标1的速度方差','目标2的速度方差',4)
end
