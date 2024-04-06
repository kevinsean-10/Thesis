function z = sobolgen(N,nVar,LB,UB)
%% generate sobol sequence in [-10,10]x[-10,10]

%%%Implementation of different sampling methods. You can program any arbitrary sampling method in this part%%%%%
%%% Here is Sobol sequences method in MATLAB %%%
p = sobolset(nVar);%,'Skip',1e4,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
A=net(p,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% The generated numbers are transferred to their correct range of parameters and rounded to two decimal places %%%
for i=1:nVar
    z(:,i)=round((LB(i)+(UB(i)-LB(i)).*A(:,i))*100)/100;
end