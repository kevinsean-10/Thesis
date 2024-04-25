function [XC,xc_clust,cl_near,XC_rad,xc_rad] = clustering(y,xcc)

[cl_num,nVar]=size(xcc);

%Find the shortest distance
for i=1:cl_num
    d(i)=norm((y-xcc(i,:)),2);
end
[Minval,cl_near]=min(d);
XC=xcc(cl_near,:);
% cl_near

% set mid-point
Xt=0.5*(y+XC);
% y
% % 
% A=FObj(Xt)
% B=FObj(y)
% C=FObj(XC)

xc_clust=[];
xc_rad=0;

if (FObj(Xt)< FObj(y)) && (FObj(Xt) < FObj(XC))
    xc_clust=y;
    xc_rad=norm((y-Xt),2);
elseif (FObj(Xt)> FObj(y)) && (FObj(Xt) > FObj(XC))
    xc_clust=y;
    xc_rad=norm((y-Xt),2);
    [XC,xc_clust,cl_near]=clustering(Xt,xcc);
else FObj(y)>FObj(XC);
    XC=y;
end

XC_rad=norm((y-Xt),2);


    
    
    


