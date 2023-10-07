close all
st=1;
num=50;
theta=rad2deg(atan2(v_y_ds,v_x_ds));
theta=theta+180;
span=0:st:360;
hist=zeros(length(span),1);
sr_dir=zeros(length(span),num);
[temp,idx]=sort(sum(spike_rate(:,2:end)),'descend');
serial=idx(1:num);
for i=1:length(theta)
    hist(floor(theta(i)/st)+1)=hist(floor(theta(i)/st)+1)+1;
    for j=1:num
        sr_dir(floor(theta(i)/st)+1,j)=sr_dir(floor(theta(i)/st)+1,j)+spike_rate(i,j+1);
    end
end
sr_dir=sr_dir(1:end-1,:);
hist=hist(1:end-1,:);
for i=1:num
   sr_dir(:,i)=sr_dir(:,i)./hist ;
end

%%
c0=ones(num,3);
for j=1:num
    for i=1:100
        c = lsqcurvefit ('Cosfit', c0(j,:), span(1:end-1),sr_dir(:,j)) ;
        c0(j,:) = c; %以计算出的 c为初值进行迭代;
    end


end

theta_preference=zeros(num,1);
for i=1:num
    theta_preference(i)=-1*c0(i,2);
    if c0(i,1)<0
        theta_preference(i)=theta_preference(i)+180;
    end
    if theta_preference(i)<0
        theta_preference(i)=theta_preference(i)+360;
    end
    
end
%%
fig=[1,4,5,7,8,9,10,15,16,20];
figure()
counter=0;
for i=1:num
    
    if ismember(i,fig)
       counter=counter+1;
       subplot(2,5,counter)
       bar(span(1:end-1),sr_dir(:,i));
       hold on
       plot(span(1:end-1),c0(i,1).*cosd(span(1:end-1)+c0(i,2))+c0(i,3));
       title('\theta =',num2str(theta_preference(i)))
       ylabel('Spike\_rate')
       xlabel('Degree')
    end
end

%%
mama=zeros(50,2);
for i=1:num
    mama(i,:)=[cosd(theta_preference(i)),sind(theta_preference(i))];
    
end
figure()
scatter(mama(:,1),mama(:,2));
xticks(-1:1)
yticks(-1:1)
xlabel('cos\theta')
ylabel('sin\theta')
% title('20160921')


