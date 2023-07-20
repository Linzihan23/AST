function ssl = SSI(X,Y)
% u1=mean2(X);
% u1=u1*u1;
% std1=std2(X);
% var1=std1*std1;
% x1=u1/var1;
% 
% u2=mean2(Y);
% u2=u2*u2;
% std3=std2(Y);
% var2=std3*std3;
% y1=u2/var2;

u1=mean2(X);
std1=std2(X);
x1=u1/std1;

u2=mean2(Y);
std3=std2(Y);
y1=u2/std3;

ssl=x1/y1;
end