clc; clear all; %close all;

%% generate sobol sequence in [-10,10]x[-10,10]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
N=250;           % number of sample points
nVar=2;           % number of decision variables


%%% UpperBand and lowerBand of each parameter%%%%
LB=[-10,-10];  UB=[10,10];
z = sobolgen(N,nVar,LB,UB);
 
%%%Implementation of different sampling methods. You can program any arbitrary sampling method in this part%%%%%
% %%% Here is Sobol sequences method in MATLAB %%%
% p = sobolset(nVar,'Skip',1e4,'Leap',1e2);
% p = scramble(p,'MatousekAffineOwen');
% A=net(p,N);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% %%% The generated numbers are transferred to their correct range of parameters and rounded to two decimal places %%%
% x1=round((LB1+(UB1-LB1).*A(:,1))*100)/100;
% x2=round((LB2+(UB2-LB2).*A(:,2))*100)/100;
% 
% plot(z(:,1),z(:,2),'.');
% pause
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fungsi Objektif masalah memaksimumkan
for i=1:N
%     z(i,:)=[x1(i),x2(i)];
    FO(i)=FObj(z(i,:));
end

%% clustering dan spiral parameter
mcl=N;
gamma=0.2;
eps=1e-7;
delta=0.01;
kcl=10;
m=250;
r=0.95;
theta=pi/4;
kmax=250;
k=1;
[Maxval,ind]=max(FO);

cl_num=1;
xc(cl_num,1:nVar)=[z(ind,1),z(ind,2)];
cl_rad(cl_num)=min(0.5*abs(UB(1)-LB(1)),0.5*(UB(2)-LB(2)));

%%
for k=1:kcl
    for i=1:mcl
         
        % Cek titik pusat cluster dan memenuhi ambang batas
        clust_center=0;
        for j=1:cl_num
            if norm((z(i,:)-xc(j,:)),2)< delta
                clust_center=1;
            end
        end
        if FObj(z(i,:))>gamma && clust_center == 0
            [XC,xc_clust,cl_near,XC_rad,xc_rad] = clustering(z(i,:),xc);
            if ~isempty(xc_clust)
                cl_num=cl_num+1;
                xc(cl_num,1:nVar)=xc_clust;
                cl_rad(cl_num)=xc_rad;
            end
%             if ~isequaln(XC,xc(cl_near,1:nVar))
                xc(cl_near,1:nVar)=XC;
                cl_rad(cl_near)=XC_rad;
%             end
        
%         i
%         [xc cl_rad']
        end
        
%          plot(z(:,1),z(:,2),'.');
%          for j=1:cl_num
%             hold on
%             plot(z(i,1),z(i,2),'o');
%             plot(xc(j,1),xc(j,2),'*');
%             viscircles(xc(j,:),cl_rad(j));
%             hold off 
%         end
%         pause(0.5)
    end

    plot(z(:,1),z(:,2),'.');
    for i=1:cl_num
        hold on
        plot(xc(i,1),xc(i,2),'*');
        viscircles(xc(i,:),cl_rad(i));
        hold off 
    end
%     pause

    for i=1:N
        FO(i)=FObj(z(i,:));
    end

   [Maxval,ind]=max(FO);  
    xp(k,:)=[z(ind,1),z(ind,2)];
    
    for i=1:N
        zb(i,:)=[r,0;0,r]*[cos(theta),-sin(theta);sin(theta),cos(theta)]*(z(i,:)'-xp(k,:)')+xp(k,:)';...
%             zb(i,:)=[r,0;0,r]*[cos(theta),-sin(theta);sin(theta),cos(theta)]*z(i,:)'-([r,0;0,r]*[cos(theta),-sin(theta);sin(theta),cos(theta)]-[1,0;0,1])*xp(k,:)';
    end
    
    z=zb;


end
% return

%% Mencari akar di setiap cluster
for n=1:cl_num
    N = m;           % number of sample points
    nVar=2;           % number of decision variables


    %%% UpperBand and lowerBand of each parameter%%%%
    LB=[xc(n,1)-cl_rad(n),xc(n,2)-cl_rad(n)];  UB=[xc(n,1)+cl_rad(n),xc(n,2)+cl_rad(n)];
    s = sobolgen(N,nVar,LB,UB);

    
  
    for k=1:kmax
        for i=1:N
            FOS(i)=FObj(s(i,:));
        end

       [Maxval,ind]=max(FOS);  
       xs(k,:)=[s(ind,1),s(ind,2)];

        for i=1:N
            zb(i,:)=[r,0;0,r]*[cos(theta),-sin(theta);sin(theta),cos(theta)]*(s(i,:)'-xs(k,:)')+xs(k,:)';...
    %             zb(i,:)=[r,0;0,r]*[cos(theta),-sin(theta);sin(theta),cos(theta)]*z(i,:)'-([r,0;0,r]*[cos(theta),-sin(theta);sin(theta),cos(theta)]-[1,0;0,1])*xp(k,:)';
        end
        s=zb;
    end

    rc(n,1:nVar)=xs(k,:); 
end

%% Cek akar di setiap cluster

root_count=1;
for n=1:cl_num
    rc(n,3)=FObj(rc(n,:));
    cek_F(n)=1-FObj(rc(n,:));
     if cek_F(n)< eps && (abs(rc(n,1))<=10 && abs(rc(n,2))<=10)
        akar(root_count,:)=rc(n,:)
        root_count=root_count+1;
     end
end

% return
%% Cek jarak antar akar
for i=1:root_count-1
    for j=i:root_count-1
        dist_root(i,j)=norm((akar(i,:)-akar(j,:)),2);
    end
end
akar_num=1;
for i=1:root_count-1
    cek=find(dist_root(:,i)<delta);
    if cek >= i
        akar_fin(akar_num,:)=akar(i,:);
        akar_num=akar_num+1;
    end   
end
akar_fin   





 %%
% clc;clear;
% ccc = [1,2,3;3,4,5;5,6,7;7,8,9]
% rrr = [1,2,3,4]
% 
% for i = 1:4
%     sibo = [ccc(i,:)'-rrr(i),ccc(i,:)'+rrr(i)]
%     LB=[ccc(i,1)-rrr(i),ccc(i,2)-rrr(i)]
%     UB=[ccc(i,1)+rrr(i),ccc(i,2)+rrr(i)]
% end