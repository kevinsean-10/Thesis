clc; clear all; %close all;

%% generate sobol sequence in [-10,10]x[-10,10]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
N=10;           % number of sample points
nVar=2;           % number of decision variables


cen=[-6.03498032104750,-0.169212173955818,0.368709371310691;
    -3.46801540386618,-0.287310386291416,0.128101196447594;
    0.0371383825640539,11.0742680147646,0.340038782800984;
    1.26036962587233,0.866932519177265,0.189580584072092;
    0.328953536131191,6.12118883403277,0.128585377804275;
    -1.75224892276408,0.394539669834870,0.193487595247569;
    0.324496999816637,-0.612637537253634,0.519200825520982;
    -0.301206037393371,4.69216346576380,0.107777161960930;
    -11.1992231276790,-0.00875241106148739,0.969795581858637;
    -4.56432333724570,0.274213189341998,0.255205915984860;
    0.110114630005849,6.56300908199895,0.510036143333613;
    -0.0699064640211145,3.49406067660502,0.301269765770983;
    -13.1387027600646,-0.0295523613451772,0.969795581858637;
    -0.976811604032224,-0.951334297159502,0.345920241436574;
    -0.477679536802524,0.590627267675109,0.533482951157123;
    -0.0320259053474273,8.65695799701224,0.547709944822630;
    -0.250396184760426,5.28336000689720,0.193088129836816;
    -0.308331156849588,4.47672693396349,0.107777161960930]

h=2;

%%% UpperBand and lowerBand of each parameter%%%%
LB=[cen(h,1)-cen(h,3),cen(h,2)-cen(h,3)];  UB=[cen(h,1)+cen(h,3),cen(h,2)+cen(h,3)]; 
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
eps=10^-7;
delta=0.01;
kcl=10;
m=250;
r=0.95;
theta=pi/4;
kmax=250;
k=1;
[Maxval,ind]=max(FO);

for k=1:kcl

    for i=1:N
        FO(i)=FObj(z(i,:));
    end

   [Maxval,ind]=max(FO);  
    xp(k,:)=[z(ind,1),z(ind,2)];
    
    for i=1:N
        zb(i,:)=[r,0;0,r]*[cos(theta),-sin(theta);sin(theta),cos(theta)]*(z(i,:)'-xp(k,:)')+xp(k,:)';
%             zb(i,:)=[r,0;0,r]*[cos(theta),-sin(theta);sin(theta),cos(theta)]*z(i,:)'-([r,0;0,r]*[cos(theta),-sin(theta);sin(theta),cos(theta)]-[1,0;0,1])*xp(k,:)';
    end
    

    plot(z(:,1),z(:,2),'.r','MarkerSize',10);
    axis([LB(1) UB(1) LB(2) UB(2)])
    grid
    hold on
    plot(zb(:,1),zb(:,2),'.b','MarkerSize',10);
%     plot(xp(k,1),xp(k,2),'og','MarkerSize',10);
    viscircles(cen(h,1:2),cen(h,3));
    hold off
    pause
    z=zb;

end
