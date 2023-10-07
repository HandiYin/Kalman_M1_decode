close all;

[le,~]=size(spike_rate);
num_test=floor(le/10);
num_train=le-num_test;
num_test=300;

X=[pos_x_ds(1:num_train)-mean(pos_x_ds(1:num_train)),pos_y_ds(1:num_train)-mean(pos_y_ds(1:num_train)),v_x_ds(1:num_train)-mean(v_x_ds(1:num_train)),v_y_ds(1:num_train)-mean(v_y_ds(1:num_train)),a_x_ds(1:num_train)-mean(a_x_ds(1:num_train)),a_y_ds(1:num_train)-mean(a_y_ds(1:num_train))]';

Z=spike_rate(1:num_train,2:end)';
X_1=X(:,1:end-1);
X_2=X(:,2:end);
A=X_2*X_1'*(X_1*X_1')^-1;
H=Z*X'*(X*X')^-1;

W=(X_2-A*X_1)*(X_2-A*X_1)'/(num_train-1);
Q=(Z-H*X)*(Z-H*X)'/num_train;


Z=spike_rate(1:num_train+num_test,2:end)';
P=eye(6);
P_minus=eye(6);
x_minus=zeros(6,num_train+num_test);
x_fix=zeros(6,num_train+num_test);
x_fix(:,1)=X(:,1);
for i=2:(num_train+num_test)
    x_minus(:,i)=A*x_fix(:,i-1);
    P_minus=A*P*A'+W;
    K=P_minus*H'*(H*P_minus*H'+Q)^-1;
    x_fix(:,i)=x_minus(:,i)+K*(Z(:,i)-H*x_minus(:,i));
    P=(eye(6)-K*H)*P_minus;    
end
%%

SNR=zeros(6,0);
figure()
subplot(2,3,1)
plot(pos_x_ds(num_train+1:num_train+num_test),'b')
hold on
plot(x_fix(1,num_train+1:num_train+num_test)+mean(pos_x_ds(1:num_train)),'r')
legend('gt','predict')
title('p_x')
R=R_square((x_fix(1,num_train+1:num_train+num_test)+mean(pos_x_ds(1:num_train))),pos_x_ds(num_train+1:num_train+num_test)');
SNR(1)=-10*log10(1-R);

subplot(2,3,4)
plot(pos_y_ds(num_train+1:num_train+num_test),'b')
hold on
plot(x_fix(2,num_train+1:num_train+num_test)+mean(pos_y_ds(1:num_train)),'r')
legend('gt','predict')
title('p_y')
R=R_square((x_fix(2,num_train+1:num_train+num_test)+mean(pos_y_ds(1:num_train))),pos_y_ds(num_train+1:num_train+num_test)');
SNR(2)=-10*log10(1-R);

subplot(2,3,2)
plot(v_x_ds(num_train+1:num_train+num_test),'b')
hold on
plot(x_fix(3,num_train+1:num_train+num_test)+mean(v_x_ds(1:num_train)),'r')
legend('gt','predict')
title('v_x')
R=R_square(x_fix(3,num_train+1:num_train+num_test)+mean(v_x_ds(1:num_train)),v_x_ds(num_train+1:num_train+num_test)');
SNR(3)=-10*log10(1-R);

subplot(2,3,5)
plot(v_y_ds(num_train+1:num_train+num_test),'b')
hold on
plot(x_fix(4,num_train+1:num_train+num_test)+mean(v_y_ds(1:num_train)),'r')
legend('gt','predict')
title('v_y')
R=R_square(x_fix(4,num_train+1:num_train+num_test)+mean(v_y_ds(1:num_train)),v_y_ds(num_train+1:num_train+num_test)');
SNR(4)=-10*log10(1-R);

subplot(2,3,3)
plot(a_x_ds(num_train+1:num_train+num_test),'b')
hold on
plot(x_fix(5,num_train+1:num_train+num_test)+mean(a_x_ds(1:num_train)),'r')
legend('gt','predict')
title('a_x')
R=R_square(x_fix(5,num_train+1:num_train+num_test)+mean(a_x_ds(1:num_train)),a_x_ds(num_train+1:num_train+num_test)');
SNR(5)=-10*log10(1-R);


subplot(2,3,6)
plot(a_y_ds(num_train+1:num_train+num_test),'b')
hold on
plot(x_fix(6,num_train+1:num_train+num_test)+mean(a_y_ds(1:num_train)),'r')
legend('gt','predict')
title('a_y')
R=R_square(x_fix(6,num_train+1:num_train+num_test)+mean(a_y_ds(1:num_train)),a_y_ds(num_train+1:num_train+num_test)');
SNR(6)=-10*log10(1-R);
