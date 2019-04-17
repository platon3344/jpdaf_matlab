%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ������  :�������ϸ������ݹ����㷨ʵ�����������˶�Ŀ��ĵ㼣�뺽���Ĺ���
% �������  ��
%           target_position:  Ŀ��ĳ�ʼλ��
%           Pd:  ������
%           Pg�� ��ȷ��������������ڵĸ��� 
%           g_sigma�� ����
%           lambda�� �Ӳ��ܶ�
%           target_delta=100�� �۲�������׼��
%           process_delta=2��  ����������׼��
%           T_target�� ��ʵĿ������Ϊ2ʱ��Ч��
%           - n��  �������� 
%           - T�� �������
%           -MC_number���������
% �������  ��
%           ��           
% ����     ��zx
% ����     �� 05-20-2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function JPDAF(target_position,Pd,Pg,g_sigma,lambda,target_delta,process_delta,T_target,n,T,MC_number)
tic

A = [1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];                   %״̬ת�ƾ���
H = [1 0 0 0;0 0 1 0];                                   %�۲����  
R=[target_delta^2  0;0  target_delta^2];                 %�۲�Э�������
Q=[process_delta^2 0;0 process_delta^2];                    %��������Э����
G=[T^2/2 0;T 0;0 T^2/2;0 T];                             %������������
P(:,:,1)=[target_delta^2 target_delta^2 0 0;target_delta^2 2*target_delta^2 0 0;0 0 target_delta^2 target_delta^2;0 0 target_delta^2 2*target_delta^2];
P(:,:,2)=[target_delta^2 target_delta^2 0 0;target_delta^2 2*target_delta^2 0 0;0 0 target_delta^2 target_delta^2;0 0 target_delta^2 2*target_delta^2];
                                                         %��ʼЭ�������

x_filter(:,:,1)=zeros(4,n);      %�洢��һ��Ŀ��ĸ�ʱ�̵��˲�ֵ
x_filter(:,:,2)=zeros(4,n);      %�洢�ڶ���Ŀ��ĸ�ʱ�̵��˲�ֵ
x_filter1(:,:,1)=zeros(4,n);     %MC_number��montle carlo��������ȫ������ۼӴ洢����
x_filter1(:,:,2)=zeros(4,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����Ŀ���ʵ��λ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
data_measurement1=zeros(4,n,2);              	%data_measurement1ʵ��λ�þ���
data_measurement1(:,1,1)=target_position(1,:)';
data_measurement1(:,1,2)=target_position(2,:)';
for tt=1:T_target
    for i=2:n
        data_measurement1(:,i,tt)=A*data_measurement1(:,i-1,tt)+G*sqrt(Q)*(randn(2,1));   %ʵ��λ�� 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%����·��
data_measurement=zeros(2,n,T_target);                  %�۲�洢����,n��������
for tt=1:T_target
    for i=1:n
       data_measurement(:,i,tt)=H*data_measurement1(:,i,tt)+sqrt(R)*randn(2,1); 
%data_measurement�۲�ֵ����,�������۲��λ��	
    end
end
for M=1:MC_number
NOISE=[];
% �˲���ʼ
for t=1:n
    y1=[];
    for tt=1:T_target
        Noise=[];
        if t~=1
            x_predic(:,tt) = A*x_filter(:,t-1,tt); %��ǰһʱ�̵��˲�ֵ��Ԥ�⵱ǰ��ֵ 
        else
            x_predic(:,tt) = target_position(tt,:)'; %��һ�β�������ʵλ�õ�Ԥ��ֵ 
        end
        P_predic(:,:,tt) = A*P(:,:,tt)*A'+G*Q*G';
        Z_predic(:,tt) = H*x_predic(:,tt);
        S(:,:,tt) = H*P_predic(:,:,tt)*H'+ R;
        K(:,:,tt) = P_predic(:,:,tt)*H'*inv(S(:,:,tt));                         %����
        ellipse_Volume(tt)=pi*g_sigma*sqrt(det(S(:,:,tt))); %�������������������������   
        number_returns(tt)=floor(10*ellipse_Volume(tt)*lambda+1); %����ز���
        side(tt)=sqrt(10*ellipse_Volume(tt))/2; %��������б߳��Ķ���֮һ
        Noise(1,:,tt)=data_measurement(1,t,tt)+side(tt)-2*rand(1,number_returns(tt))*side(tt); 
													%��Ԥ��ֵ��Χ��������ز�
        Noise(2,:,tt)=data_measurement(2,t,tt)+side(tt)-2*rand(1,number_returns(tt))*side(tt);
        NOISE=[NOISE Noise(:,:,tt)];
        %
        b=zeros(2,1);
        b=data_measurement(:,t,tt);
        y1=[y1 Noise(:,:,tt) b];            %�����յ������еĻز�����y1��
    end
    y=[];
    d=[];
    Ohm=[];                                            %ȷ�Ͼ���
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
           m=m+1;                        %��¼�ܵ���Ч�ز�����
           y=[y y1(:,j)];              %����������������е����лز�����y��
        end
    end
    
    %m=0��ʾ����Ч�ز�
         
    if m==0                                %�޻ز������
       x_filter(:,t,1)= x_predic(:,1);
       x_filter(:,t,2)=x_predic(:,2);
       P(:,:,1)=P_predic(:,:,1);                                           
       P(:,:,2)=P_predic(:,:,2);
    else
       if m*2-sum(Ohm(:))==0              %��ͬĿ�����Ч�ز�û�й��������
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
                    Bk=lambda*sqrt(2*pi*det(S(:,:,tt)))*(1-Pd*Pg)/Pd;      %��b0 
                    E=zeros(1,m);
                    belta=zeros(1,m);
                    for i=1:m
                        a=(yy(:,i)-Z_predic(:,tt))'*inv(S(:,:,tt))*(yy(:,i)-Z_predic(:,tt));
                        E(i)=E(i)+exp(-a/2);
                    end 
                    belta0=Bk/(Bk+sum(E));             %�޻ز�ʱ�Ĺ�������
                    v=zeros(2,1);
                    v1=zeros(2,2);
                    for i=1:m
                        belta(i)=E(i)/(Bk+sum(E));                 %���������
                        v=v+belta(i)*(yy(:,i)-Z_predic(:,tt));
                        v1=v1+belta(i)*(yy(:,i)-Z_predic(:,tt))*(yy(:,i)-Z_predic(:,tt))';
                    end
                    x_filter(:,t,tt)= x_predic(:,tt) + K(:,:,tt)*v;
                    %��Э����
                    Pc=(eye(4)-K(:,:,tt)*H)*P_predic(:,:,tt);
                    PP=K(:,:,tt)*(v1-v*v')*K(:,:,tt)';
                    P(:,:,tt)=belta0*P_predic(:,:,tt)+(1-belta0)*Pc+PP;
               end
           end
       else
                                             		%��ͬĿ���лز����������
           [Ohm_predict,Nry]=AResolution(Ohm);        %��ȷ�Ͼ����ֳɻ�������
           c=0;                                     %��һ��ϵ��
           beta=zeros(m,T_target);                      %��������
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
                %��Э����
                Pc=(eye(4)-K(:,:,tt)*H)*P_predic(:,:,tt);
                PP=K(:,:,tt)*(v1-v*v')*K(:,:,tt)';
                P(:,:,tt)=beta0*P_predic(:,:,tt)+(1-beta0)*Pc+PP;
           end
       end
    end
end
x_filter1(:,:,1)=x_filter1(:,:,1)+x_filter(:,:,1);       %ÿ�����ؿ���ʵ��õ����˲�ֵ���ۼ�
x_filter1(:,:,2)=x_filter1(:,:,2)+x_filter(:,:,2);
end
toc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %��ͼ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_filter(:,:,1)=x_filter1(:,:,1)/MC_number;                         %�˲�ֵ��ƽ��
x_filter(:,:,2)=x_filter1(:,:,2)/MC_number;
figure                               %����Ŀ��Ĺ���λ�ã�����λ�ã��Ӳ�λ��
color=['r-' 'r:' 'k-' 'k:'];
for tt=1:T_target
plot(x_filter(1,:,tt),x_filter(3,:,tt),color(2*tt-1:2*tt)),hold on
%plot(data_measurement1(1,:),data_measurement1(3,:),'-')
plot(data_measurement1(1,:,tt),data_measurement1(3,:,tt),color(2*tt+3:2*tt+4))
%axis([0 30 1 25])
end
plot(NOISE(1,:),NOISE(2,:),'.')                               %�Ӳ�λ��
xlabel('x(m)'),ylabel('y(m)');
legend('Ŀ��1�Ĺ���λ��','Ŀ��1����ʵλ��','Ŀ��2�Ĺ���λ��','Ŀ��2����ʵλ��','�Ӳ���λ��',4)

for tt=1:T_target
    figure                                %����Ŀ��X��Ĺ���λ�ã�����λ��
    plot(1:n,x_filter(1,:,tt),'r-'),hold on
    plot(data_measurement1(1,:,tt),'k-')
    title(sprintf('Ŀ��%d',tt));
    xlabel('t(s)'),ylabel('x(m)');
    legend('X�������λ��','X������ʵλ��',4)
end

for tt=1:T_target
    figure                                %����Ŀ��Y��Ĺ���λ�ã�����λ��
    plot(1:n,x_filter(3,:,tt),'r-'),hold on
    plot(data_measurement1(3,:,tt),'k-')
    title(sprintf('Ŀ��%d',tt));
    xlabel('t(s)'),ylabel('y(m)');
    legend('Y�������λ��','Y������ʵλ��',4)
end

figure											%����Ŀ��λ�ù��ƾ������
color=['r:' 'b-'];
for tt=1:T_target
    a=zeros(1,n);                                   %λ�þ������ 
    for j=1:n
a(1,j)=sqrt((x_filter(1,j,tt)-data_measurement1(1,j,tt))^2+(x_filter(3,j,tt)-data_measurement1(3,j,tt))^2);
    end
    plot(1:n,a(1,:),color(2*tt-1:2*tt));hold on;
    xlabel('t(s)'),ylabel('Ԥ�����(m)');
end
legend('Ŀ��1�����','Ŀ��2�����',4)

figure											%����Ŀ��λ�ù��ƾ������
color=['r:' 'b-'];
for tt=1:T_target  
    a=zeros(1,n);                                    %�ٶȾ������
    for j=1:n
a(1,j)=sqrt((x_filter(2,j,tt)-data_measurement1(2,j,tt))^2+(x_filter(4,j,tt)-data_measurement1(4,j,tt))^2);
    end
    plot(1:n,a(1,:),color(2*tt-1:2*tt));hold on;
    xlabel('t(s)'),ylabel('Ԥ�����(m/s)');
end
    legend('Ŀ��1���ٶȷ���','Ŀ��2���ٶȷ���',4)
end
