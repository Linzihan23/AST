function nc = NC(X,Y)
   m1=mean2(X);
   m2=mean2(Y);
   a=X-m1;
   b=Y-m2;
   c=(a.*b);
   d=a.*a;
   e=b.*b;
   f=sum(c(:));
   g=sqrt((sum(d(:))*sum(e(:))));
   nc=f/g;
  
end