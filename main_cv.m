clc
close all

% global fu fv cu cv h c1 c2 s1 s2

image = im2double(imread('f00378.png'));
image = rgb2gray(image);
[mo,no] = size(image);
% fu = 700.1; fv = 700.1;
fu = 35; fv = 35*9/16;
cu = 669.3; cv = 339.0;
h = 1706;
pitch = 1.0; yaw = 0.0;
c1 = cos(pitch); c2 = cos(yaw);
s1 = sin(pitch); s2 = sin(yaw);

% Icropped = image(400:end,:);

% IPM = ipm_function(image);

% image_hor = image(339:end,:);

figure(1)
imshow(image)
% hold on
% plot( 1280, 720, '.g', 'LineWidth', 3);
% plot( 1280, 1, '.g', 'LineWidth', 3);
% 
% 
% 
% for i  = 1:1280
% image(339,i) = 0;
% end
% 
% figure(2)
% imshow(image)
% 
% figure(3)
% imshow(image_hor)

T = h*[(-1/fu)*c2 (1/fv)*s1*s2 ((1/fu)*cu*c2)-((1/fv)*cv*s1*s2)-(c1*s2) 0;... 
    (1/fu)*s2 (1/fv)*s1*c1 (-(1/fu)*cu*s2)-((1/fv)*cv*s1*c2)-(c1*c2) 0;...
    0 (1/fv)*c1 ((-1/fv)*cv*c1)+s1 0;...
    0 (-1/(h*fv))*c1 ((1/(h*fv))*cv*c1)-s1/h 0];

image_hor = image(400:end,:);
[m,n] = size(image_hor);
proj = zeros(m,n,2);
for i=400:mo
    for j = 1:no
        comp_proj = T*[j;i;1;1];
        comp_proj = comp_proj/(comp_proj(4));
        proj(i-399,j,1) = comp_proj(1);
        proj(i-399,j,2) = comp_proj(2);
%         figure 
%         hold on
%         plot(proj(i,j,1),proj(i,j,2))
    end
end 
% figure
% % proj1 = 
% scatter(proj(:,:,1),proj(:,:,2))
ipm = zeros(240,360);
[m1,n1] = size(ipm);
k = 0;
proj_reshape = reshape(proj,321*no,2);
xmax = max(proj_reshape(:,1));
xmin = min(proj_reshape(:,1));
ymax = max(proj_reshape(:,2));
ymin = min(proj_reshape(:,2));

for i=400:mo
    for j = 1:no
        k = k+1;
        ipmx = ((n1-1)/(xmax-xmin))*(proj(i-399,j,1)-xmin)+1;
        ipmy = ((m1-1)/(ymin-ymax))*(proj(i-399,j,2)-ymax)+1;
        ipmx = round(ipmx); ipmy = round(ipmy);
        ipm(ipmy,ipmx) = image(i,j);
    end
end
Xvis = proj(:,:,1); Yvis = proj(:,:,2);
figure(2)
imshow(ipm)
% plot(Xvis(:), Yvis(:), 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 4);

% yg = (0.1:0.001:0.66);
% xg = (-1.95:0.01:1.77)';
% 
% 
% Iw = griddata(Xvis, Yvis, Icropped, xg, yg, 'linear');
% figure;
% imagesc(Iw,[0 1]); 
% colormap(gray);
% axis image;

%% Detecting vertical edges
% im_filtered = imfilter(ipm, [-1 ;0 ;1]');
im_filtered = imgaussfilt(ipm,'FilterSize',[3,3]);
[Gx,Gy] = gradient(im_filtered);

[im_filtered,Gxy] = gradient(Gx);
figure
imshow(Gxy)

%% Thresholding
im_thres = im_filtered;
gray_mean = mean(im_thres);
gray_mean = 0.975*mean(gray_mean);
indices = find(abs(im_filtered)<gray_mean);
im_thres(indices) = 0;
% cutoff_frequency = 1;
% filter = fspecial('Gaussian', cutoff_frequency*4+1, cutoff_frequency);
% low_freq = imfilter(im_thres,filter);
figure
imshow(im_thres)