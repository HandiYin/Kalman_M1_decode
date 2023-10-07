%% 载入数据
clear all;
close all;
load('indy_20160921_01.mat');

%% 计算spike rate，筛选神经细胞
[m,n]=size(spikes);
timetable=t(1):0.064:t(end); %64毫秒时间窗
total_t=t(end)-t(1);
timetable=timetable';
A=zeros(length(timetable),1);
B=[];
for i=1:m
    for j=2:n
        temp=spikes{i,j};
        if ~isempty(temp)
            if length(temp)/total_t >0.5 %平均发放率大于0.5Hz
                for k=1:length(timetable)-1
                    temp_freq=length(temp(temp>timetable(k)&temp<=timetable(k+1)));
                    A(k)=temp_freq;
                end
                B=[B,A];
            end
           
        end
    end  
end

spike_rate=B/0.064;

spike_rate=[timetable,spike_rate];
save indy_20160921_01_SR.mat spike_rate

%% 计算运动学参数
fs=250;%Hz
finger_2D=finger_pos*[0,0;-10,0;0,-10];
fin_pos_x=finger_2D(:,1);
fin_pos_y=finger_2D(:,2);
[B,A]=butter(4,10*2/250,'low');
fin_pos_x_filter=filter(B,A,fin_pos_x);
fin_pos_y_filter=filter(B,A,fin_pos_y);
% figure()
% subplot(2,1,1)
% plot(fin_pos_x)
% subplot(2,1,2)
% plot(fin_pos_x_filter)
fin_v_x=diff(fin_pos_x_filter)*fs;% mm/s
fin_v_y=diff(fin_pos_y_filter)*fs;
fin_a_x=diff(fin_v_x)*fs;% mm/s^2
fin_a_y=diff(fin_v_y)*fs;
pos_x_ds=resample(fin_pos_x_filter,1000,250*64);
pos_y_ds=resample(fin_pos_y_filter,1000,250*64);
v_x_ds=resample(fin_v_x,1000,250*64);
v_y_ds=resample(fin_v_y,1000,250*64);
a_x_ds=resample(fin_a_x,1000,250*64);
a_y_ds=resample(fin_a_y,1000,250*64);
% figure()
% subplot(2,1,1)
% plot(fin_pos_x)
% subplot(2,1,2)
% plot(pos_x_ds)

%% 线性回归 start！
% 留下后1/10点作为测试集，前9/10的点作为训练集
close all;
[le,~]=size(spike_rate);
num_test=floor(le/10);
num_train=le-num_test;
num_test=300;
train_p_x=[ones(num_train,1),spike_rate(1:num_train,2:end)];
[b_p_x,bint_p_x,r_p_x,rint_p_x,stats_p_x] = regress(pos_x_ds(1:num_train),train_p_x);
test_p_x=[ones(num_test,1),spike_rate(num_train+1:num_train+num_test,2:end)]*b_p_x;
figure()
subplot(2,3,1)
plot(pos_x_ds(num_train+1:num_train+num_test),'b')
hold on
plot(test_p_x,'r')
legend('gt','predict')
title('p_x')


train_v_x=[ones(num_train,1),spike_rate(1:num_train,2:end)];
[b_v_x,bint_v_x,r_v_x,rint_v_x,stats_v_x] = regress(v_x_ds(1:num_train),train_v_x);
test_v_x=[ones(num_test,1),spike_rate(num_train+1:num_train+num_test,2:end)]*b_v_x;
subplot(2,3,4)
plot(v_x_ds(num_train+1:num_train+num_test),'b')
hold on
plot(test_v_x,'r')
legend('gt','predict')
title('v_x')


train_a_x=[ones(num_train,1),spike_rate(1:num_train,2:end)];
[b_a_x,bint_a_x,r_a_x,rint_a_x,stats_a_x] = regress(a_x_ds(1:num_train),train_a_x);
test_a_x=[ones(num_test,1),spike_rate(num_train+1:num_train+num_test,2:end)]*b_a_x;
subplot(2,3,2)
plot(a_x_ds(num_train+1:num_train+num_test),'b')
hold on
plot(test_a_x,'r')
legend('gt','predict')
title('a_x')

train_p_y=[ones(num_train,1),spike_rate(1:num_train,2:end)];
[b_p_y,bint_p_y,r_p_y,rint_p_y,stats_p_y] = regress(pos_y_ds(1:num_train),train_p_y);
test_p_y=[ones(num_test,1),spike_rate(num_train+1:num_train+num_test,2:end)]*b_p_y;
subplot(2,3,5)
plot(pos_y_ds(num_train+1:num_train+num_test),'b')
hold on
plot(test_p_y,'r')
legend('gt','predict')
title('p_y')

train_v_y=[ones(num_train,1),spike_rate(1:num_train,2:end)];
[b_v_y,bint_v_y,r_v_y,rint_v_y,stats_v_y] = regress(v_y_ds(1:num_train),train_v_y);
test_v_y=[ones(num_test,1),spike_rate(num_train+1:num_train+num_test,2:end)]*b_v_y;
subplot(2,3,3)
plot(v_y_ds(num_train+1:num_train+num_test),'b')
hold on
plot(test_v_y,'r')
legend('gt','predict')
title('v_y')

train_a_y=[ones(num_train,1),spike_rate(1:num_train,2:end)];
[b_a_y,bint_a_y,r_a_y,rint_a_y,stats_a_y] = regress(a_y_ds(1:num_train),train_a_y);
test_a_y=[ones(num_test,1),spike_rate(num_train+1:num_train+num_test,2:end)]*b_a_y;
subplot(2,3,6)
plot(a_y_ds(num_train+1:num_train+num_test),'b')
hold on
plot(test_a_y,'r')
legend('gt','predict')
title('a_y')




