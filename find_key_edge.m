clear
useresult=1;  %是否用预测的8_BA500_4t002nL1B1degreepredictionedge8_SW500_8_05t01nL1B1degreepredictionedge
if useresult
%     y=load('8_SW500_8_05t01nL1B1degreepredictionedge.mat','maxatkline');
    y=load('8_BA500_4t002nL1B1degreepredictionedge.mat','maxatkline');
%     y=load('6_BA1000_3t002nL1B1degreepredictionedge.mat','maxatkline');
    
    maxatkline=y.maxatkline;
attack=1;
destn=1;
n1=1;
repair=0;
refine=0;
Lfindex=0;
t=0.08;
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
Be=textread('F:\py\Mtb\8_BA500_4\edge_betweenness_centrality2file.txt','', 'headerlines', 0);% DU

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
% for i=1:1:m
%         if Ss(L(i,2))~=0
%             exAS(L(i,1),L(i,2))=(A(L(i,1),L(i,2))/(Ss(L(i,2))));
%             
%         end
% end
for i=1:1:m
        if Ss(L(i,2))~=0
%             exAS(L(i,1),L(i,2))=max(A(L(i,1),L(i,2))/(Ss(L(i,2))),A(L(i,1),L(i,2))/(Ss(L(i,1))));
% exAS(L(i,1),L(i,2))=(A(L(i,1),L(i,2))/(Ss(L(i,2))));
           exAS(L(i,1),L(i,2))=(A(L(i,1),L(i,2))/(Ss(L(i,2)))+A(L(i,1),L(i,2))/(Ss(L(i,1))));
%             exAS(L(i,1),L(i,2))=(A(L(i,1),L(i,2))/(S(L(i,2)))+A(L(i,1),L(i,2))/(S(L(i,1))));
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
[resum2,rerows2]=ismember(cascLfail(:,1:2),Be(:,[1 2]),'rows');

cascLfail=[cascLfail,maxatkline(rerows,3),exASorder(rerows3,3),Be(rerows2,3)]; 
for ii=1:length(cascLfail(:,1))
loada=Loadc(cascLfail(ii,1),cascLfail(ii,2));
cascLfail(ii,8)=loada;
end
avgcascLfail=mean(cascLfail,2);
cascLfailpercent=size(find(cascLfail(:,4)==1),1)/cascnumber;
cascLfailpercentavg=mean(cascLfail(:,4));
end
%save('cascattack8_BA1000_4t01nL02B1degree.mat','cascLfail','cascLfailpercent','LFaillinecell','Lfailcell','erfacell','cascavgerfa','avgcascLfail');
end
% [resum1,rerows1]=ismember(cascLfail(find(cascLfail(:,4)==1),1:2),exASorder(:,[1 2]),'rows');
% [resum,rerows]=ismember(cascLfail(find(cascLfail(:,4)==1),1:2),maxatkline(:,[1 2]),'rows');
cascLfailtheta=cascLfail(find(cascLfail(:,4)==1),:);     
% plot(maxatkline(rerows,3),cascLfail(find(cascLfail(:,4)==1),3),'m.')
% title('4_SW100_4_05t005')
fprintf('结束')
% plot(cascLfail(:,5),cascLfail(:,4))
save('ADD_Edge_betweenness_pic8_BA500_4t008nL1B1degreepredictionedgeEsasort.mat','cascLfailtheta','cascLfail');
% save('ADD_Edge_betweenness_pic6_BA1000_3t010nL1B1degreepredictionedgeEsasort.mat','cascLfailtheta','cascLfail');
% save('Edge_betweenness_pic8_SW500_8_05t008nL1B1degreepredictionedgeEsasort.mat','cascLfailtheta','cascLfail');
% save('pic4_SW100_4_08t02nL1B1degreepredictionedge.mat','cascLfailtheta','cascLfail');
index=1;
cx=cascLfail(:,4);%result
cy=cascLfail(:,5);%Triggered the first stage neighbouring edge maximized 
cz=cascLfail(:,6);%Edge weight ratio
cbe=cascLfail(:,7).^index
cload=cascLfail(:,8)

figure(1)
subplot(2,1,1)
plot(cy,cx,'kd')
% xlabel('Triggered the first stage neighbouring edge maximized α')

% xlabel('${\alpha}_{\Gamma_{ij}}^{max}$ of a removal edge $l_{ij}$','Interpreter','latex')
xlabel({'The threshold value ${\alpha}_{\Gamma_{ij}}^{max}$ of a removal edge $l_{ij}$'},'Interpreter','latex')



ylabel('The avalanche size $P_f $','Interpreter','latex')
title(['Spearman correlation = ',num2str(corr(cy,cx, 'type' , 'Spearman'))])
subplot(2,1,2)
plot(cz,cx,'kd')
xlabel('Edge weight ratio $e_{ij}$ of a removal edge $l_{ij}$ ','Interpreter','latex')
ylabel('The avalanche size $P_f $','Interpreter','latex')
legend('The avalanche size $P_f $','Interpreter','latex')
title(['Spearman correlation = ',num2str(corr(cz,cx, 'type' , 'Spearman'))])
figure(3)
% subplot(2,1,1)
 plot(cz,cy,'b.','LineWidth',2)
 xlabel('Edge weight ratio $e_{ij}$ of a removal edge $l_{ij}$ ','Interpreter','latex')
 
%  xlabel('The value $w_{ij}k_{j}/s_{j}+w_{ij}k_{i}/s_{i} $ of a removal edge $l_{ij}$ ','Interpreter','latex')

ylabel({'The threshold value ${\alpha}_{\Gamma_{ij}}^{max}$ of a removal edge $l_{ij}$'},'Interpreter','latex')

legend('${\alpha}_{\Gamma_{ij}}^{max}$','Interpreter','latex')
title(['Spearman correlation = ',num2str(corr(cz,cy, 'type' , 'Spearman'))])
%  plot(cbe,cz,'kd')
% subplot(2,1,2)Pearson


figure(4)
cbe=cascLfail(:,7);
plot(cbe,cz,'b.');
% end
xlabel('Edge betweenness centrality $ EBC_{ij} $ of a removal edge $l_{ij}$','Interpreter','latex')
ylabel('Edge weight ratio $e_{ij}$ of a removal edge $l_{ij}$','Interpreter','latex')
% legend('The avalanche size $P_f $')
% title(['Spearman correlation = ',num2str(corr(cbe,cz, 'type' , 'Spearman'))])
title(['Pearson correlation = ',num2str(corr(cbe,cz, 'type' , 'Pearson'))])

% plot(cascLfailtheta(:,4),cascLfailtheta(:,5))
figure(5)
cx1=(cascLfail(find(cascLfail(:,4)==1),6));
limy=cz
zb = (min(limy):(max(limy)-min(limy))/50:max(limy));
[countscx,centerscx] = hist(cx1, zb);
[countscz,centerscz] = hist(limy, zb);
bar(centerscx, countscx ./ countscz,'FaceColor','b')
xlabel('Edge weight ratio $e_{ij}$ of a removal edge $l_{ij}$','Interpreter','latex')
ylabel('The complete breakdown ratio ${P_r}$','Interpreter','latex')
legend('Probability')

figure(6)
% plot(cy,cz)
cx2=(cascLfail(find(cascLfail(:,4)==1),5));
limx=cy
zb1 = (min(limx):(max(limx)-min(limx))/50:max(limx));
[countscx1,centerscx1] = hist(cx2, zb1);
[countscz1,centerscz1] = hist(limx, zb1);
bar(centerscx1, countscx1 ./ countscz1,'FaceColor','b')
% xlabel('${\alpha}_{\Gamma_{ij}}^{max}$ of a removal edge $l_{ij}$','Interpreter','latex')
xlabel({'The threshold value ${\alpha}_{\Gamma_{ij}}^{max}$ of a removal edge $l_{ij}$'},'Interpreter','latex')

ylabel('The complete breakdown ratio ${P_r}$','Interpreter','latex')
legend('Probability')

figure(7)
% plot(cy,cz)
cx3=(cascLfail(find(cascLfail(:,4)==1),7));
limbe=cbe
zb2 = (min(limbe):(max(limbe)-min(limbe))/50:max(limbe));
[countscx1,centerscx1] = hist(cx3, zb2);
[countscz1,centerscz1] = hist(limbe, zb2);
bar(centerscx1, countscx1 ./ countscz1,'FaceColor','b')
xlabel('Edge betweenness centrality $ EBC_{ij} $ of a removal edge $l_{ij}$','Interpreter','latex')
ylabel('The complete breakdown ratio ${P_r}$','Interpreter','latex')
legend('Probability')

figure(8)
% plot(cy,cz)
cx4=(cascLfail(find(cascLfail(:,4)==1),8));
limload=cload
zb3 = (min(limload):(max(limload)-min(limload))/50:max(limload));
[countscx1,centerscx1] = hist(cx4, zb3);
[countscz1,centerscz1] = hist(limload, zb3);
% subplot(2,1,1) ,[.69 .69 .69]
bar(centerscx1, countscx1 ./ countscz1,'FaceColor','b')
xlabel('Initial load of a removal edge $l_{ij}$','Interpreter','latex')
ylabel('The complete breakdown ratio ${P_r}$','Interpreter','latex')
legend('Probability')
figure(10)
bar(centerscz1,countscz1)
figure(9)
subplot(2,1,1)
plot(cbe,cx,'kd')
xlabel('Edge betweenness centrality $ EBC_{ij} $ of a removal edge $l_{ij}$','Interpreter','latex')
ylabel('The avalanche size $P_{f}$','Interpreter','latex')
legend('The avalanche size $P_f $','Interpreter','latex')
title(['Spearman correlation = ', num2str(corr(cbe,cx, 'type' , 'Spearman'))])
subplot(2,1,2)
plot(cload,cx,'kd')
xlabel('Initial load of a removal edge $l_{ij}$','Interpreter','latex')
ylabel('The avalanche size $P_{f}$','Interpreter','latex')
legend('The avalanche size $P_f $','Interpreter','latex')
title(['Spearman correlation = ',num2str(corr(cload,cx, 'type' , 'Spearman'))])
figure(11)
% subplot(2,1,1)
 plot(cbe,cload,'b.','LineWidth',2)
 xlabel('Edge weight ratio $e_{ij}$ of a removal edge $l_{ij}$ ','Interpreter','latex')
%  xlabel('The value $w_{ij}k_{j}/s_{j}+w_{ij}k_{i}/s_{i} $ of a removal edge $l_{ij}$ ','Interpreter','latex')
ylabel({'Initial load of a removal edge $l_{ij}$'},'Interpreter','latex')

title(['Spearman correlation = ',num2str(corr(cbe,cload, 'type' , 'Spearman'))])