function enl = ENL(img_roi)
    % ENL Summary of this function goes here 
    % img_roi is the selected regions from the your image(or the image )
    % img_roi=double(img_roi);
    
u1=mean2(img_roi);
u1=u1*u1;
std1=std2(img_roi);
var1=std1*std1;
enl=u1/var1;
end