clear
load Data190
x2=Data;
load GroundTruth(10class)
x1=GroundTruth;
load Training_Samples_Image_1
x=Training_Samples_Image;
ca=11;band=190;s=2;
save 'd2.mat' x2
[a,b,c1]=size(x2);
x4=x1(1:144,1:144)+1;%data=x2;
data(:,:,1:190)=x2(:,:,1:190);
%   data(:,:,11:20)=x2(:,:,50:59);
%  data(:,:,31:60)=x2(:,:,101:130);
nx=b-1;                 %Number of colloms
ny=a-1;                  %Number of lines
nb=band;                  %Number of spectral bands
nc=ca;                   %Number of clesses
nt=50;
ni(1:17)=0;
train=zeros(nb,nt,nc);%pq(145,145)=0;
for k=1:16
for iy=1:145
    for ix=1:145
        if x(iy,ix)==k
            ni(k)=ni(k)+1;train3(:,ni(k),k)=data(iy,ix,:);%pq(iy,ix)=k;
        end
    end
end
end

for iy=1:6:ny
for ix=1:3:nx
if x1(iy,ix)==0 && ni(17)<nt;
	        ni(17)=ni(17)+1;
            train3(:,ni(17),17)=data(iy,ix,:);%pq(iy,ix)=17;
end
end
end
% 1,4,7,9,13,16  1,7,9
train(:,:,1)=train3(:,:,17);
i=0;train(:,:,2-i)=train3(:,:,2);train(:,:,3-i)=train3(:,:,3);
train(:,:,4-i)=train3(:,:,5);train(:,:,5-i)=train3(:,:,6);train(:,:,6-i)=train3(:,:,8);
train(:,:,7-i)=train3(:,:,10);train(:,:,8-i)=train3(:,:,11);train(:,:,9-i)=train3(:,:,12);
train(:,:,10-i)=train3(:,:,14);train(:,:,11-i)=train3(:,:,15);ni=[50,50,50,50,50,50,50,50,50,50,50];
save 'train.mat' train
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meanc=zeros(nb,nc);  %mean vector of all classes
for ic=1:nc
    meanc(:,ic)=(mean(train(:,1:ni(ic),ic)'))';
end   

% % cxi=zeros(nb,nb,nc);  %inverse of covariance matrix of classes
% % deti=zeros(nc,1);     %log of determinan of covariance matrix
% % for ic=1:nc
% % cxi(:,:,ic)=(train(:,1:ni(ic),ic)*train(:,1:ni(ic),ic)')/(ni(ic)-1)-meanc(:,ic)*meanc(:,ic)';
% %     deti(ic)=det(cxi(:,:,ic));
% %     cxi(:,:,ic)=inv(cxi(:,:,ic));
% %     deti(ic)=log(deti(ic));
% % end
% % % caculation the measure
dist=zeros(nc,1);
xi=zeros(nb,1);
yi=zeros(nb,1);
thh=zeros(nc,1);
for ic=1:nc
   for it=1:ni(ic)
         xi(:)=train(:,it,ic)-meanc(:,ic);
        % th(ic,it)=deti(ic)+xi'*cxi(:,:,ic)*xi;
         th(ic,it)=norm(xi,2);
  end
end
thh=max(th');           %thershold of class rejection
  gt_map=x4;
class_map=zeros(ny,nx);
for iy=1:ny
   for ix=1:nx
      xi(:)=data(iy,ix,:);
      for ic=1:nc
% %           yi=xi-meanc(:,ic);
% %           dist(ic)=deti(ic)+yi'*cxi(:,:,ic)*yi;
      dist(ic)=norm(xi-meanc(:,ic),2);
      end
      [ma nu]=min(dist);
      if ma<=thh(nu)
          class_map(iy,ix)=nu;
      else
         class_map(iy,ix)=1;
      end
      end
end


     figure; subplot(1,4,1); imagesc(x1);xlabel('gtmap')
             subplot(1,4,2); imagesc(class_map);xlabel('classmap')
             
             con=zeros(nc+1,nc+3);       %make the confusion matrix
for iy=1:ny
   for ix=1:nx
      i=gt_map(iy,ix);
      j=class_map(iy,ix);
      con(i,j)=con(i,j)+1;
   end
end
con(nc+1,:)=sum(con);
con(:,nc+1)=(sum(con'))';
%calculation of accuracy
con(nc+1,nc+2)=0;
con(nc+1,nc+1)=0;
con(nc+1,nc+3)=0;
for ic=1:nc
   con(ic,nc+2)=(con(ic,ic)/con(ic,nc+1))*100;
   con(ic,nc+3)=(con(ic,ic)/con(nc+1,ic))*100;
   con(nc+1,nc+2)=con(nc+1,nc+2)+con(ic,ic);
   con(nc+1,nc+1)=con(nc+1,nc+1)+con(ic,nc+1);
   con(nc+1,nc+3)=con(nc+1,nc+3)+con(nc+1,ic);
end
con(nc+1,nc+3)=(100*con(nc+1,nc+2)/con(nc+1,nc+3));
con(nc+1,nc+2)=(100*con(nc+1,nc+2)/con(nc+1,nc+1));
con1=double(con);

[a,b]=size(class_map);conn=double(con)/(a*b);
j=1:b;
for i=1:a
        y1((i-1)*b+j)=gt_map(i,j);
end
t=1:nc;
for i=1:nc
    z1(i)=conn(i,i);
end
z2=histc(y1,t);z2=z2/size(y1,2);
z=[z1;z2];
zm=min(z);
kappa=(sum(zm)-sum(z1.*z2))/(1-sum(z1.*z2));
%**** in this section calculate confusion matrix in the edge pixels for class_map
x5=class_map;
x4=gt_map;
e1=0;cone=zeros(nc+1,nc+3);n=0;y4(a,b)=0;
for i=1:a
    for j=1:b
    i1=s*fix((i+s-1)/s)-(s-1);j1=s*fix((j+s-1)/s)-(s-1);
if ne(sum(sum(abs(x4(i1:i1+s-1,j1:j1+s-1)/x4(i1,j1)-ones(s,s)))),0);
    n=n+1;cone(x4(i,j),x5(i,j))=cone(x4(i,j),x5(i,j))+1;y3(n)=x4(i,j);
if ne(x4(i,j),x5(i,j))
    e1=e1+1;
end
    end
    end
end
cone(nc+1,:)=sum(cone);
cone(:,nc+1)=(sum(cone'))';
%calculation of accuracy in the edge pixels of class_map
cone(nc+1,nc+2)=0;
cone(nc+1,nc+1)=0;
cone(nc+1,nc+3)=0;
for ic=1:nc
   cone(ic,nc+2)=(cone(ic,ic)/cone(ic,nc+1))*100;
   cone(ic,nc+3)=(cone(ic,ic)/cone(nc+1,ic))*100;
   cone(nc+1,nc+2)=cone(nc+1,nc+2)+cone(ic,ic);
   cone(nc+1,nc+1)=cone(nc+1,nc+1)+cone(ic,nc+1);
   cone(nc+1,nc+3)=cone(nc+1,nc+3)+cone(nc+1,ic);
end
cone(nc+1,nc+3)=(100*cone(nc+1,nc+2)/cone(nc+1,nc+3));
cone(nc+1,nc+2)=(100*cone(nc+1,nc+2)/cone(nc+1,nc+1));
e1=1-e1/n;con1e=cone;conen=cone/n;
t=1:nc;
for i=1:nc
    z1(i)=conen(i,i);
end
z2=histc(y3,t);z2=z2/size(y3,2);
z=[z1;z2];
zm=min(z);
kappa1=(sum(zm)-sum(z1.*z2))/(1-sum(z1.*z2));
%**************************************************
[con(12,14)/100 ; kappa;e1; kappa1] %pcc & kappa for class_map

 load d2h;
s2=double(x);s2=imresize(s2,s,'nearest');s2=(s2-1.6)*255/241.4;
data=s2;
class_map1=zeros(ny,nx);
for iy=1:ny
   for ix=1:nx
      xi(:)=data(iy,ix,:);
      for ic=1:nc
          yi=xi-meanc(:,ic);
          dist(ic)=deti(ic)+yi'*cxi(:,:,ic)*yi;
      end
      [ma nu]=min(dist);
      if ma<=thh(nu)
         class_map1(iy,ix)=nu;
      else
         class_map1(iy,ix)=nc;
      end
      end
end

con=zeros(nc+1,nc+3);       %make the confusion matrix for hard classification
for iy=1:ny
   for ix=1:nx
      i=gt_map(iy,ix);
      j=class_map1(iy,ix);
      con(i,j)=con(i,j)+1;
   end
end
con(nc+1,:)=sum(con);
con(:,nc+1)=(sum(con'))';
%calculation of accuracy for hard classification
con(nc+1,nc+2)=0;
con(nc+1,nc+1)=0;
con(nc+1,nc+3)=0;
for ic=1:nc
   con(ic,nc+2)=(con(ic,ic)/con(ic,nc+1))*100;
   con(ic,nc+3)=(con(ic,ic)/con(nc+1,ic))*100;
   con(nc+1,nc+2)=con(nc+1,nc+2)+con(ic,ic);
   con(nc+1,nc+1)=con(nc+1,nc+1)+con(ic,nc+1);
   con(nc+1,nc+3)=con(nc+1,nc+3)+con(nc+1,ic);
end
con(nc+1,nc+3)=(100*con(nc+1,nc+2)/con(nc+1,nc+3));
con(nc+1,nc+2)=(100*con(nc+1,nc+2)/con(nc+1,nc+1));
[a,b]=size(class_map1);conn=double(con)/(a*b);conh=con;
j=1:b;
for i=1:a
        y1((i-1)*b+j)=gt_map(i,j);
end
t=1:nc;
for i=1:nc
    z1(i)=conn(i,i);
end
z2=histc(y1,t);z2=z2/size(y1,2);
z=[z1;z2];
zm=min(z);
kappa=(sum(zm)-sum(z1.*z2))/(1-sum(z1.*z2));
%**in this section calculate confusion matrix in the edge pixels for the
%hardclassification
x5=class_map1;
x4=gt_map;
e1=0;cone=zeros(nc+1,nc+3);n=0;y4(a,b)=0;
for i=1:a
    for j=1:b
    i1=s*fix((i+s-1)/s)-(s-1);j1=s*fix((j+s-1)/s)-(s-1);
if ne(sum(sum(abs(x4(i1:i1+s-1,j1:j1+s-1)/x4(i1,j1)-ones(s,s)))),0);
    n=n+1;cone(x4(i,j),x5(i,j))=cone(x4(i,j),x5(i,j))+1;y3(n)=x4(i,j);
if ne(x4(i,j),x5(i,j))
    e1=e1+1;
end
    end
    end
end
cone(nc+1,:)=sum(cone);
cone(:,nc+1)=(sum(cone'))';
%calculation of accuracy for edge pixle of hardclassification
cone(nc+1,nc+2)=0;
cone(nc+1,nc+1)=0;
cone(nc+1,nc+3)=0;
for ic=1:nc
   cone(ic,nc+2)=(cone(ic,ic)/cone(ic,nc+1))*100;
   cone(ic,nc+3)=(cone(ic,ic)/cone(nc+1,ic))*100;
   cone(nc+1,nc+2)=cone(nc+1,nc+2)+cone(ic,ic);
   cone(nc+1,nc+1)=cone(nc+1,nc+1)+cone(ic,nc+1);
   cone(nc+1,nc+3)=cone(nc+1,nc+3)+cone(nc+1,ic);
end
cone(nc+1,nc+3)=(100*cone(nc+1,nc+2)/cone(nc+1,nc+3));
cone(nc+1,nc+2)=(100*cone(nc+1,nc+2)/cone(nc+1,nc+1));
e1=1-e1/n;conhe=cone;conen=cone/n;
t=1:nc;
for i=1:nc
    z1(i)=conen(i,i);
end
z2=histc(y3,t);z2=z2/size(y3,2);
z=[z1;z2];
zm=min(z);
kappa1=(sum(zm)-sum(z1.*z2))/(1-sum(z1.*z2));
%**************************************************
[con(18,20)/100 ; kappa;e1; kappa1]   %for hard classification
subplot(1,4,3); imagesc(class_map1);xlabel('classmap1')

load d2h4;
s2=double(x);s2=(s2-1.6)*255/241.4;
data=s2;
class_map2=zeros(ny,nx);
for iy=1:ny
   for ix=1:nx
      xi(:)=data(iy,ix,:);
      for ic=1:nc
          yi=xi-meanc(:,ic);
          dist(ic)=deti(ic)+yi'*cxi(:,:,ic)*yi;
      end
      [ma nu]=min(dist);
      if ma<=thh(nu)
         class_map2(iy,ix)=nu;
      else
         class_map2(iy,ix)=nc;
      end
      end
end
             subplot(1,4,4); imagesc(class_map2);xlabel('classmap2')

con=zeros(nc+1,nc+3);       %make the confusion matrix for proposed method
for iy=1:ny
   for ix=1:nx
      i=gt_map(iy,ix);
      j=class_map2(iy,ix);
      con(i,j)=con(i,j)+1;
   end
end
con(nc+1,:)=sum(con);
con(:,nc+1)=(sum(con'))';
%calculation of accuracy for proposed method
con(nc+1,nc+2)=0;
con(nc+1,nc+1)=0;
con(nc+1,nc+3)=0;
for ic=1:nc
   con(ic,nc+2)=(con(ic,ic)/con(ic,nc+1))*100;
   con(ic,nc+3)=(con(ic,ic)/con(nc+1,ic))*100;
   con(nc+1,nc+2)=con(nc+1,nc+2)+con(ic,ic);
   con(nc+1,nc+1)=con(nc+1,nc+1)+con(ic,nc+1);
   con(nc+1,nc+3)=con(nc+1,nc+3)+con(nc+1,ic);
end
con(nc+1,nc+3)=(100*con(nc+1,nc+2)/con(nc+1,nc+3));
con(nc+1,nc+2)=(100*con(nc+1,nc+2)/con(nc+1,nc+1));
conp=double(con);
[a,b]=size(class_map2);conn=double(con)/(a*b);
j=1:b;
for i=1:a
        y1((i-1)*b+j)=gt_map(i,j);
end
t=1:nc;
for i=1:nc
    z1(i)=conn(i,i);
end
z2=histc(y1,t);z2=z2/size(y1,2);
z=[z1;z2];
zm=min(z);
kappa=(sum(zm)-sum(z1.*z2))/(1-sum(z1.*z2));
%*%**in this section calculate confusion matrix in the edge pixels for the
%proposed method
x5=class_map2;
x4=gt_map;
e1=0;cone=zeros(nc+1,nc+3);n=0;y4(a,b)=0;
for i=1:a
    for j=1:b
    i1=s*fix((i+s-1)/s)-(s-1);j1=s*fix((j+s-1)/s)-(s-1);
if ne(sum(sum(abs(x4(i1:i1+s-1,j1:j1+s-1)/x4(i1,j1)-ones(s,s)))),0);
    n=n+1;cone(x4(i,j),x5(i,j))=cone(x4(i,j),x5(i,j))+1;y3(n)=x4(i,j);
if ne(x4(i,j),x5(i,j))
    e1=e1+1;
end
    end
    end
end
cone(nc+1,:)=sum(cone);
cone(:,nc+1)=(sum(cone'))';
%calculation of accuracy for edge pixels of proposed method
cone(nc+1,nc+2)=0;
cone(nc+1,nc+1)=0;
cone(nc+1,nc+3)=0;
for ic=1:nc
   cone(ic,nc+2)=(cone(ic,ic)/cone(ic,nc+1))*100;
   cone(ic,nc+3)=(cone(ic,ic)/cone(nc+1,ic))*100;
   cone(nc+1,nc+2)=cone(nc+1,nc+2)+cone(ic,ic);
   cone(nc+1,nc+1)=cone(nc+1,nc+1)+cone(ic,nc+1);
   cone(nc+1,nc+3)=cone(nc+1,nc+3)+cone(nc+1,ic);
end
cone(nc+1,nc+3)=(100*cone(nc+1,nc+2)/cone(nc+1,nc+3));
cone(nc+1,nc+2)=(100*cone(nc+1,nc+2)/cone(nc+1,nc+1));
e1=1-e1/n;conpe=cone;conen=cone/n;
t=1:nc;
for i=1:nc
    z1(i)=conen(i,i);
end
z2=histc(y3,t);z2=z2/size(y3,2);
z=[z1;z2];
zm=min(z);
kappa1=(sum(zm)-sum(z1.*z2))/(1-sum(z1.*z2));
%**************************************************
[con(18,20)/100 ; kappa;e1; kappa1] %for proposed method