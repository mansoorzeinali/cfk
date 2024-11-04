close all hidden
clear
x=importdata('datanew2.mat');  x3=importdata('gt_map2.mat');ca2=3;ca1=3;s1=12;s2=8;th1=.9992;n1=90;n2=100;%fuzzy,correlation
tic
p(1)=fix(size(x,1)./s1).*s1;p(2)=fix(size(x,2)./s2).*s2;
a1=p(1);b1=p(2);c1=size(x,3);
a2=a1/s1;b2=b1/s2;k3=0;
x5=double(x(1:a1,1:b1,1:c1));
x6=double(x3(1:a1,1:b1,:));
for k1=1:s1
    for k2=1:s2
    g=x5((k1-1)*a2+1:k1*a2,(k2-1)*b2+1:k2*b2,:);k3=k3+1;
    j=1:b2;
for i=1:a2
        yz((i-1)*b2+j,1:c1)=g(i,j,1:c1);
end
 opts = [nan;200;nan;0];
[c5,h1,u1]=fcm(yz,ca1,opts);c6((k3-1)*ca1+1:k3*ca1,:)=c5;%s4=h1;h1=h1';h2(:,:,k3)=h1;
    end
end
t=a1:-1:1;
t1=1:b1;
figure,pcolor(t1,t,x6(1:a1,1:b1));axis equal tight
shading flat
       [h3,c7]=kmeans(c6,ca2,'distance','cosine','MaxIter',1000,'Replicates',5);qa=-1;h3=h3';
  c7(ca2+1,1:c1)=zeros(1,c1);
    x8(a1,b1)=1;
for i=1:a1
    for j=1:b1
        w6=000000000;
        for l=1:ca2
            w5=permute(x5(i,j,:),[3,2,1]);
          w7=((c7(l,:)*w5))/((norm(w5')*norm(c7(l,:))));
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
