clear
useresult=1;  %是否用预测的8_BA500_4t002nL1B1degreepredictionedge8_SW500_8_05t01nL1B1degreepredictionedge
if useresult
%     y=load('8_SW500_8_05t01nL1B1degreepredictionedge.mat','maxatkline');
    y=load(' 8_BA500_4t002nL1B1degreepredictionedge.mat','maxatkline');
    maxatkline=y.maxatkline;
attack=1;
destn=1;
n1=1;
repair=0;
refine=0;
Lfindex=0;
t=0.05;
cycle=1; %0==求每个边的theta,1==级联故障。8_BA500_4
cascattack=1;  %级联故障
cascatkrefine=0;
attackcycle=0; %攻击找边界最大值
reeverytheta=[];
% B=textread('F:\py\pyscrptis\nx.betweenness_centrality.txt','' , 'headerlines', 0);% 中介度
B=textread('F:\py\Mtb\8_BA500_4\nx.degree_centrality.txt','' , 'headerlines', 0);% 中介度
% B=textread('F:\py\pyscrptis\nx.degree_centrality.txt','' , 'headerlines', 0);% 中介度
% B=textread('F:\py\pyscrptis\nx.closeness_centrality.txt','' , 'headerlines', 0);% 中介度
B(find(B~=0))=1e+5*(B(find(B~=0))).^(1);
% D=textread('F:\py\pyscrptis\node_degree.txt','', 'headerlines', 0);% DU
D=textread('F:\py\Mtb\8_BA500_4\node_degree.txt','', 'headerlines', 0);% DU
v=size(B,1);
A=zeros(v,v);
Load=zeros(v,v);
Z1=ones(v,v);
% L=textread('F:\py\pyscrptis\edgelist.txt','', 'headerlines', 0);%%连接关系
L=textread('F:\py\Mtb\8_BA500_4\edgelist.txt','', 'headerlines', 0);%%连接关系
% Be=textread('F:\py\Mtb\8_BA500_4\edge_betweenness_centrality2file.txt','', 'headerlines', 0);% DU

% 寻找连接边是否存在
[m,n]=size(L);
for i=1:1:m
    A(L(i,1),L(i,2))=[D(L(i,1))*D(L(i,2))]^n1;
end
A=A+A';
%M=mean(A,2);%每一行平均值
S=sum(A,2);
Ss=S./sum(A~=0,2);
% Ss=S;
[explorevalue exploresort]=sort(B./S);
BS=B./S;
for i=1:1:m
    BSline(i,1)=BS(L(i,1));
    BSline(i,2)=BS(L(i,2));
end

exBSorder=[explorevalue exploresort];
exBSReverseorder=flipud(exBSorder);
exAS=[];
for i=1:1:m
        if Ss(L(i,2))~=0
            exAS(L(i,1),L(i,2))=(A(L(i,1),L(i,2))/(Ss(L(i,2))));
        end
end

% for i=1:1:v
%     for j=1:1:v
%         if Ss(j)~=0
%             exAS(i,j)=(A(i,j)/(Ss(j)));
%         end
%     end
% end
[exASi,exASj]=find(exAS~=0);
exASorder=sortrows([exASi,exASj,exAS(exAS~=0)],-3);

for j=1:1:v
    if S(j)~=0
        for i=1:1:v
            Load(j,i)=A(j,i)/S(j)*B(j);
        end
    end
end
Loadc=Load;
F2=zeros(v,v);
for j=1:1:v
    for i=1:1:v
        if (B(i)*S(j)+B(j)*S(i))~=0
            F2(j,i)=S(i)*S(j)/(B(i)*S(j)+B(j)*S(i));
        end
    end
end
if cascattack
    erfacell={};
    Lfailcell={};
    LFaillinecell={};
    cascavgerfa=[];
    cascLfail=[];
    cascnumber=0;
    cascam=0;
    Lcascatktimes=[];
    for cascatk=1:1:size(maxatkline,1)
        erfa=[];
        A2=A;
%          exasorder=0;
%         cascatki=maxatkline(cascatk,1);
%         cascatkj=maxatkline(cascatk,2);
% %         exasorder=1;
        cascatki=exASorder(cascatk,1);
        cascatkj=exASorder(cascatk,2);
        
        Lcascatk=[cascatki,cascatkj];
        A2(cascatki,cascatkj)=0;
        A2(cascatkj,cascatki)=0;
        if cascatkrefine&&ismember(Lcascatk,reeverythetarow(:,1:2),'rows')
            cascam=cascam+1;
            [cascatkjnumber,cascatkjrows]=ismember(Lcascatk,reeverythetarow(:,1:2),'rows');
            cascatkires=reeverythetarow(cascatkjrows,3);
            cascatkjres=reeverythetarow(cascatkjrows,4);
            res=(D(cascatkires)*D(cascatkjres))^n1;
            A2(cascatkires,cascatkjres)=res;
            A2(cascatkjres,cascatkires)=res;
%           修复两边
            cascatkires=reeverythetarow(cascatkjrows-1,3);
            cascatkjres=reeverythetarow(cascatkjrows-1,4);
            res=(D(cascatkires)*D(cascatkjres))^n1;
            A2(cascatkires,cascatkjres)=res;
            A2(cascatkjres,cascatkires)=res;
            fprintf('循环修复')
        end
        Lfail=1;
        Lfailline=[];
        cyclenumber=0;
        Lfaillinecell={};
        cascnumber=cascnumber+1;
        S=sum(A2,2);
        while cycle&any(any(A2))
            S=sum(A2,2);
            F0=A2;
            n=0;
            cyclenumber=cyclenumber+1;
            F1=zeros(v,v);
            F0((F0~=0))=1;
            for j=1:1:v
                for i=1:1:v
                    if (S(j)*S(i))~=0
                        F1(j,i)=F0(j,i)*[Z1(j,i)*F2(j,i)*[B(j)/S(j)+B(i)/S(i)]-1];
                        if (F1(j,i)>t);
                            A2(j,i)=0;
                            n=n+1;
                        end
                    end
                end
            end
            erfa=[erfa,max(max(triu(F1)))];
%     [Lfaili,Lfailj]=find(triu(F1)>0);
% Lfaillinecell{cyclenumber}=sortrows([Lcascatk,Lfaili,Lfailj,F1(triu(F1)>0),A(sub2ind(size(A),Lfaili,Lfailj))],-5);
%             Lfailline=[Lfailline;Lfaili,Lfailj];
            Lfail=[Lfail,n/2];
            if max(max(F1))<=t;
                break;
            end
            Lcascatktimes=[Lcascatktimes;Lcascatk,cyclenumber,n/2];
        end
% LFaillinecell{cascnumber}=Lfaillinecell;
%         erfacell{cascnumber}=erfa;
        Lfailcell{cascnumber}=Lfail;
%         cascavgerfa=[cascavgerfa,mean(erfa,2)];
cascLfail=[cascLfail;Lcascatk,cyclenumber,sum(Lfail,2)/m];   
    end
[resum3,rerows3]=ismember(cascLfail(:,1:2),exASorder(:,[1 2]),'rows');
[resum,rerows]=ismember(cascLfail(:,1:2),maxatkline(:,[1 2]),'rows');
cascLfail=[cascLfail,maxatkline(rerows,3),exASorder(rerows3,3)]; 
avgcascLfail=mean(cascLfail,2);
cascLfailpercent=size(find(cascLfail(:,4)==1),1)/cascnumber;
cascLfailpercentavg=mean(cascLfail(:,4));
end
%save('cascattack8_BA1000_4t01nL02B1degree.mat','cascLfail','cascLfailpercent','LFaillinecell','Lfailcell','erfacell','cascavgerfa','avgcascLfail');
end
[resum1,rerows1]=ismember(cascLfail(find(cascLfail(:,4)==1),1:2),exASorder(:,[1 2]),'rows');
[resum,rerows]=ismember(cascLfail(find(cascLfail(:,4)==1),1:2),maxatkline(:,[1 2]),'rows');
cascLfailtheta=[cascLfail(find(cascLfail(:,4)==1),1:3) maxatkline(rerows,3) exASorder(rerows1,3) ];     
% plot(maxatkline(rerows,3),cascLfail(find(cascLfail(:,4)==1),3),'m.')
% title('4_SW100_4_05t005')
fprintf('结束')
plot(cascLfail(:,5),cascLfail(:,4))
% save('pic4_BA100_2t012nL1B1degreepredictionedge.mat','cascLfailtheta','cascLfail');
% save('pic4_SW100_4_08t02nL1B1degreepredictionedge.mat','cascLfailtheta','cascLfail');
% save('pic8_BA500_4t008nL1B1degreepredictionedgeEsasort.mat','cascLfailtheta','cascLfail');%zhege
% save('pic8_BA500_4t012nL1B1degreepredictionedgeEsasort.mat','cascLfailtheta','cascLfail');%zhege
save('pic8_BA500_4t005nL1B1degreepredictionedgeEsasort.mat','cascLfailtheta','cascLfail');%zhege
% save('pic8_BA500_4t008nL1B1degreepredictionedgeEsasort.mat','cascLfailtheta','cascLfail');%zhege
% dyuiyin对应 11
% save('pic8_SW500_8_05t008nL1B1degreepredictionedgeEsasort.mat','cascLfailtheta','cascLfail');
%对应12
cx=cascLfail(:,4);
cy=cascLfail(:,5);
cz=cascLfail(:,6);
figure(1)
subplot(2,1,1)
plot(cy,cx,'kd')
xlabel('Triggered the first stage neighbouring edge maximized α')
ylabel('The avalanche size P_{f}')
title({'Spearman',num2str(corr(cy,cx, 'type' , 'Spearman'))})
subplot(2,1,2)
plot(cz,cx,'rd')
xlabel('Edge weight ratio w_{ij}k_{j}/s_{j}')
ylabel('The avalanche size P_{f}')
legend('The avalanche size $P_f $')
title({'Spearman',num2str(corr(cz,cx, 'type' , 'Spearman'))})
figure(3)
% subplot(2,1,1)
 plot(cz,cy,'b.','LineWidth',2)
 xlabel('Edge weight ratio w_{ij}k_{j}/s_{j}')
ylabel({'Triggered maximized α of ';' neighbouring edge at the first step'})
legend('Threshold value ')

% subplot(2,1,2)
% plot(cascLfailtheta(:,4),cascLfailtheta(:,5))
figure(5)
cx1=(cascLfail(find(cascLfail(:,4)==1),6));
zb = (min(cz):(max(cz)-min(cz))/50:max(cz));
[countscx,centerscx] = hist(cx1, zb);
[countscz,centerscz] = hist(cz, zb);

bar(centerscx, countscx ./ countscz,'FaceColor',[.69 .69 .69])
xlabel('Edge weight ratio w_{ij}k_{j}/s_{j}')
ylabel('Probability of edge failure')
legend('Probability')
figure(4)
% plot(cy,cz)
cx2=(cascLfail(find(cascLfail(:,4)==1),5));
zb1 = (min(cy):(max(cy)-min(cy))/50:max(cy));
[countscx1,centerscx1] = hist(cx2, zb1);
[countscz1,centerscz1] = hist(cy, zb1);

bar(centerscx1, countscx1 ./ countscz1,'FaceColor',[.69 .69 .69])
xlabel('Triggered the first stage neighbouring edge maximized α')
ylabel('Probability of edge failure')
legend('Probability')
% yy=hist(sApPro,x); %计算各个区间的个数cx1=(cascLfail(find(cascLfail(:,4)==1),6));
% yy=yy/length(sApPro); %计算各个区间的比例
% bar(x,yy);
% reshape(cascLfail(find(cascLfail(:,4)==1),3),9,3)
% nn=4;
% ynumber=[];
% if  exasorder==1
%     j=5;
% else
%     j=4;
% end
% sApProMax=max(cascLfailtheta(:,j));
% sApProMin=min(cascLfailtheta(:,j));
% x=linspace(sApProMin,sApProMax+0.1,nn+1);
% cascLfailnumberunique=sort(unique(cascLfailtheta(:,1)),1);
% for i=1:1:nn
% [fx,fy]=find(cascLfailtheta(:,j)>=x(i)&cascLfailtheta(:,j)<x(i+1));
% cascLfailnumber=cascLfailtheta(fx,1);
% yy=hist(cascLfailtheta(fx,1),cascLfailnumberunique);
% if ~any(cascLfailnumber)
% yy=yy';
% end
% ynumber=[ynumber;yy];
% end
% 
% bar3(ynumber,'grouped')
% set(gca, 'XTickLabel', cascLfailnumberunique)
% set(gca, 'YTickLabel', x)
% legend(num2str(cascLfailnumberunique))
% zlabel('Periodic numbers')
% title('The correlation between variable factor and Cascading failure of Periodic numbers')
% if  exasorder==1
% ylabel('Edge weight ratio')
% else
% ylabel('Edge first-triggered maxmum threshhold ')
% end



% h=bar3(ynumber,'grouped');
% h1=get(h,'children');
% set(h1{1},'FaceColor','r')
% text ([1 1]',[4 5-1/6]',[83 87]',num2str([83 87]'),'HorizontalAlignment','center','VerticalAlignment','bottom')
% title('Grouped Style')
% figure
% bar3(Z,'stacked')
% title('Stacked Style')
% a1=cascLfail(:,5);
% b1=(cascLfail(:,6));
% c1=(cascLfail(:,4));
% e1=hist(c1);




