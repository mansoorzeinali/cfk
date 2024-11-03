close all hidden
clear
% x10=imread('bin.bmp');x=im2double(x10(:,:,1:3));x3=ones(size(x));ca2=2;ca1=2;s=80;th1=.994;%fuzzy,correlation;n1=35;n2=45;n1=55;n2=65;ca1=2
% x=importdata('Salinas_corrected.mat');x3=importdata('Salinas_gt.mat');ca2=16;ca1=16;s1=32;s2=11;%ca1=3;s=64;th1=.9988;
% load('Data190.mat');x=Data(1:128,1:128,:); load('GroundTruth(10class).mat');x3=GroundTruth(1:128,1:128,:);ca2=10;ca1=11;s1=32;s2=32;th1=.991;%fuzzy,correlation
%  load PaviaU.mat; load PaviaU_gt.mat;x=paviaU(:,:,:);x3=paviaU_gt(:,:);ca2=9;ca1=9;s1=23;s2=17;%x=x/max(max(max(x)));
%  R: band 55, G: band 33, and B: band 13
%   load Pavia.mat; load Pavia_gt.mat;x=pavia(:,:,:);x3=pavia_gt(:,:);ca2=9;ca1=9;s1=16;s2=55;%x=x/max(max(max(x)));
% x=importdata('data12band.mat');  x3=importdata('gt12band.mat');ca2=9;ca1=9;s1=44;s2=28;th1=.9995;
% x=importdata('data.mat'); x3=importdata('gt12.mat');ca2=12;ca1=12;s1=12;s2=19;th1=.9995;
% x=importdata('SalinasA_corrected.mat'); x3=importdata('SalinasA_gt.mat');ca2=6;ca1=6;s1=8;s2=8;
% R: band 57, G: band 27, and G: band 17
% x=importdata('data.mat');  x3=importdata('gt_map4.mat');ca2=4;ca1=4;s1=20;s2=15;th1=.994;%fuzzy,correlation;n1=35;n2=45;n1=55;n2=65;ca1=2
% load('Indian_pines_corrected.mat');load('GroundTruth(10class).mat')
%  load('jasperwl');x=js1;x11=importdata('jaspermap.mat');x3=x11(2:99,2:99);s1=16;s2=16;ca2=4;th1=.996;ca1=3;
% load('jasper.mat');x=x5; x3=importdata('jaspermap.mat');ca2=4;th1=.996;ca1=4;s1=10; s2=10;%or s=16fuzzy,combinen1=40;n2=60;n1=80;n2=100;n1=120;n2=140
x=importdata('datanew2.mat');  x3=importdata('gt_map2.mat');ca2=3;ca1=3;s1=12;s2=8;th1=.9992;n1=90;n2=100;%fuzzy,correlation
% x=importdata('12003.mat');  x3=importdata('12003.mat');ca2=3;ca1=3;s1=20;s2=20;th1=.9992;n1=90;n2=100;%fuzzy,correlation
% x=imread('DC.tif');x=im2double(x);x3=importdata('ground.mat');x3=im2double(x3);ca1=8;ca2=8;s1=32;s2=20;
% data1 = imread('dc.tif');[nx, ny, nz] = size(data1);data1 = double(data1);
% img = data1 - repmat(min(min(data1,[],1),[],2),[nx ny]);img = img ./ repmat(max(max(img,[],1),[],2),[nx ny]);
% image = img;imshow(cat(3,image(:,:,60),image(:,:,27),image(:,:,17)));
% load GT.tif,[Y,newmap] = imapprox(cdata,colormap,8);% imshow(Y, newmap);
% x=imnoise(x,'gaussian',0,.001);
tic
% th1=.995;n1=35;n2=45;delta=.005;main
th1=.995;n1=75;n2=85;delta=.005;
p(1)=fix(size(x,1)./s1).*s1;p(2)=fix(size(x,2)./s2).*s2;
a1=p(1);b1=p(2);c1=size(x,3);
a2=a1/s1;b2=b1/s2;k3=0;
x5=double(x(1:a1,1:b1,1:c1));%ma=max(max(max(x5)));mi=min(min(min(x5)));x5=(x5-mi)/(ma-mi);
x6=double(x3(1:a1,1:b1,:));
%  x5=x5.*[x6~=0];
%  for i=1:a1
%      for j=1:b1
%          if x5(i,j,1)==0
%              x5(i,j,:)=ma*ones(1,c1);
%          end
%      end
%  end
% za=min(min(min(x5)));zb=max(max(max(x5)));zc=sqrt(c1)*(zb-za);
for k1=1:s1
    for k2=1:s2
    g=x5((k1-1)*a2+1:k1*a2,(k2-1)*b2+1:k2*b2,:);k3=k3+1;
    j=1:b2;
for i=1:a2
        yz((i-1)*b2+j,1:c1)=g(i,j,1:c1);
end
% m7=1;
% for    j=1:b2
% for i=1:a2
%     if sum(g(i,j,1:c1))>0
%         yz(m7,1:c1)=g(i,j,1:c1);m7=m7+1;
%     end
%  end
%  end
%    [h1,c5]=classm2(g,ca1);h2(:,:,k3)=h1;c6((k3-1)*ca1+1:k3*ca1,:)=c5;s4=h1;
%    [h1,c5]=kmeans(yz,ca1);h1=h1;h2(:,:,k3)=h1;c6((k3-1)*ca1+1:k3*ca1,:)=c5;s4=h1;
%    [h1,c5]=kmeans(yz,ca1,'distance','cosine','MaxIter',1000,'Replicates',5);h1=h1;h2(:,:,k3)=h1;c6((k3-1)*ca1+1:k3*ca1,:)=c5;s4=h1;
 opts = [nan;200;nan;0];
[c5,h1,u1]=fcm(yz,ca1,opts);c6((k3-1)*ca1+1:k3*ca1,:)=c5;%s4=h1;h1=h1';h2(:,:,k3)=h1;
%  params = struct('km_eps',.001,'km_max_iter',100);[c5,h1,u1]=it2_fcm(yz,ca1,2,2,1,.001,200,1,params);h1=h1';h2(:,:,k3)=h1;c6((k3-1)*ca1+1:k3*ca1,:)=c5;s4=h1';
% j=1:b2;
% for i=1:a2
%         x1(i,j,:)=s4((i-1)*b2+j,:);
% end
% for i1=1:a2
%     for j1=1:b2
% %           [o,x2((k1-1)*a2+i1,(k2-1)*b2+j1)]=max(x1(i1,j1,:)); %FCM
% %          x2((k1-1)*a2+i1,(k2-1)*b2+j1)=max(x1(i1,j1,:)); %KMEANS
%     end
% end
 
    end
end
t=a1:-1:1;
t1=1:b1;
%  figure,pcolor(t1,t,x2(1:a1,1:b1));axis equal tight
% shading flat
figure,pcolor(t1,t,x6(1:a1,1:b1));axis equal tight
shading flat


% % ˜ÇåÔ ÊÚÏÇÏ ãÑÇ˜Ò%
% nx=size(c6,1);ix=1;
%  while nx>n2 || nx<n1
% clear c2 w;
%     c2(1,:)=c6(1,:);l2=1;
% for i=1:size(c6,1)
%     for j=1:size(c2,1)
% %         w(j)=norm(c2(j,:)-c6(i,:));
%                 w(j)=abs(c2(j,:)*c6(i,:)')/(norm(c2(j,:))*(norm(c6(i,:))));
%     end
%                 if max(w) < th1
%             c2(l2+1,:)=c6(i,:);l2=l2+1;
%                end
% end
% nx=size(c2,1);ix=ix+1;nxi(ix)=nx;thx(ix)=th1;
% % th1=th1-.0001;
% th1=th1+((n2-nx)/abs(n2-nx))*delta;
% delta=delta/2;
% end
%        [h3,c7]=kmeans(c6,ca2,'distance','cosine','MaxIter',1000);qa=-1;h3=h3';
       [h3,c7]=kmeans(c6,ca2,'distance','cosine','MaxIter',1000,'Replicates',5);qa=-1;h3=h3';
%        [h3,c7]=kmeans(c6,ca2);qa=-1;h3=h3';
%   [c7,h3,u2]=fcm(c6,ca2,opts);qa=-1;
%  params = struct('km_eps',.001,'km_max_iter',100);[c7,h3,u2]=it2_fcm(c2,ca2,2,2,1,.001,200,1,params);
  c7(ca2+1,1:c1)=zeros(1,c1);
    x8(a1,b1)=1;
for i=1:a1
    for j=1:b1
        w6=000000000;
        for l=1:ca2
            w5=permute(x5(i,j,:),[3,2,1]);
          w7=((c7(l,:)*w5))/((norm(w5')*norm(c7(l,:))));
%             w7=norm(w5'-c7(l,:),2);
%              w8=norm(w5'-c7(l,:),2);w9=1-w8/zc;w7=mean([w9,w7]); 
            if w7>w6
                w6=w7;x8(i,j)=l;
            end
        end
    end
end
 x9=x8.*[x6~=0];
time=toc;
t=a1:-1:1;
t1=1:b1;
figure,pcolor(t1,t,x9(1:a1,1:b1));axis equal tight
shading flat
c7=c7;x1=x8;
save 'map2' a1 b1 x6 x1 ca2 time c7 x9
