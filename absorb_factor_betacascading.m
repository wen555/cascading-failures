clc;
clear;
cascLfailpercentcycleT=[];
cascLfailpercentcycleTbeta=[];
betan=0;
for beta=0:0.1:1;
    betan=betan+1;
for t=0:0.02:0.5
attack=1;
destn=1;
n1=1;
repair=0;
refine=0;
Lfindex=0;
cycle=1; %0==求每个边的theta,1==级联故障
attackandrefine=0;
recyclerefine=0;   %循环修复
cascatkrefine=0;   %级联修复
cascattack=1;  %级联故障
attackcycle=0; %攻击找边界最大值
%B=textread('F:\py\Mtb\6_BA500_3\nx.closeness_centrality.txt','' , 'headerlines', 0);% 中介度
reeverytheta=[];
%B=textread('F:\py\pyscrptis\nx.betweenness_centrality.txt','' , 'headerlines', 0);% 中介度
%B=textread('F:\py\pyscrptis\nx.degree_centrality.txt','' , 'headerlines', 0);% 中介度
B=textread('F:\py\Mtb\4_BA50_2\nx.degree_centrality.txt','' , 'headerlines', 0);% 中介度
B(find(B~=0))=1e+5*(B(find(B~=0))).^(1);
% D=textread('F:\py\pyscrptis\node_degree.txt','', 'headerlines', 0);% DU
D=textread('F:\py\Mtb\4_BA50_2\node_degree.txt','', 'headerlines', 0);% DU
v=size(B,1);
A=zeros(v,v);
Load=zeros(v,v);
Z1=ones(v,v);
% L=textread('F:\py\pyscrptis\edgelist.txt','', 'headerlines', 0);%%连接关系
L=textread('F:\py\Mtb\4_BA50_2\edgelist.txt','', 'headerlines', 0);%%连接关系
% 寻找连接边是否存在
[m,n]=size(L);
for i=1:1:m
    A(L(i,1),L(i,2))=[D(L(i,1))*D(L(i,2))]^n1;
end
A=A+A';
%M=mean(A,2);%每一行平均值
S=sum(A,2);
Ss=sum(A,2);
% [explorevalue exploresort]=sort(B./S);
% BS=B./S;
% for i=1:1:m
%     BSline(i,1)=BS(L(i,1));
%     BSline(i,2)=BS(L(i,2));
% end
% 
% exBSorder=[explorevalue exploresort];
% exBSReverseorder=flipud(exBSorder);
% for i=1:1:v
%     for j=1:1:v
%         if S(j)~=0
%             exAS(i,j)=A(i,j)/S(j);
%         end
%     end
% end
% [exASi,exASj]=find(exAS~=0);
% exASorder=sortrows([exASi,exASj,exAS(exAS~=0),A(sub2ind(size(A),exASi,exASj))],1);

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
    for cascatk=1:1:m
        erfa=[];
        A2=A;
        cascatki=L(cascatk,1);
        cascatkj=L(cascatk,2);
        Lcascatk=[cascatki,cascatkj];
        A2(cascatki,cascatkj)=0;
        A2(cascatkj,cascatki)=0;
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
                        F1(j,i)=F0(j,i)*[Z1(j,i)*F2(j,i)*beta*[B(j)/Ss(j)*(Ss(j)-S(j))/S(j)+B(i)/Ss(i)*(Ss(i)-S(i))/S(i)]];
                        if (F1(j,i)>t);
                            A2(j,i)=0;
                            n=n+1;
                        end
                    end
                end
            end
            erfa=[erfa,max(max(triu(F1)))];
            [Lfaili,Lfailj]=find(triu(F1)>t);
            Lfaillinecell{cyclenumber}=sortrows([Lfaili,Lfailj,F1(triu(F1)>t),A(sub2ind(size(A),Lfaili,Lfailj))],-3);
            Lfailline=[Lfailline;Lfaili,Lfailj];
            Lfail=[Lfail,n/2];
            if max(max(F1))<=t;
                break;
            end
        end
        LFaillinecell{cascnumber}=Lfaillinecell;
        erfacell{cascnumber}=erfa;
        Lfailcell{cascnumber}=Lfail;
        cascavgerfa=[cascavgerfa,mean(erfa,2)];
        cascLfail=[cascLfail,sum(Lfail,2)];
    end
    avgcascLfail=mean(cascLfail-1,2)/m;
    cascLfailpercent=size(find(cascLfail./m==1),2)/cascnumber;  
end
cascLfailpercentcycleT=[cascLfailpercentcycleT;t,cascLfailpercent,avgcascLfail];
end
cascLfailpercentcycleTbeta=[cascLfailpercentcycleTbeta;cascLfailpercentcycleT];
end
% cascLfailpercentcycleT=cascLfailpercentcycleT';
% save('8_BA500_4nL1B1degreenesscascadingfailT.mat','cascLfailpercentcycleT')
fprintf('over')
cascLfailpercentcycleT1=[];
i=[1:26];
figure(1)
xx=[0:0.02:0.5]';
Beta=[0:0.1:1];
for n=1:1:11
cascLfailpercentcycleT1=[cascLfailpercentcycleT1,cascLfailpercentcycleT(i+(n-1)*26*ones(1,26),3)]
end  
% plot(xx,cascLfailpercentcycleT1(:,2),xx,cascLfailpercentcycleT1(:,3),xx,cascLfailpercentcycleT1(:,4),xx,cascLfailpercentcycleT1(:,5),xx,cascLfailpercentcycleT1(:,6),xx,cascLfailpercentcycleT1(:,7),xx,cascLfailpercentcycleT1(:,8),xx,cascLfailpercentcycleT1(:,9),xx,cascLfailpercentcycleT1(:,10),xx,cascLfailpercentcycleT1(:,11));
plot(xx,cascLfailpercentcycleT1(:,2),'Color',[0   0  1],'linewidth',1.5);hold on;
plot(xx,cascLfailpercentcycleT1(:,3),'Color',[1   0  0] ,'linewidth',1.5);hold on;
plot(xx,cascLfailpercentcycleT1(:,4),'Color',[ 0   1  0 ],'linewidth',1.5);hold on;
plot(xx,cascLfailpercentcycleT1(:,5),'Color',[0   0  1 ] ,'linewidth',1.5);hold on;
plot(xx,cascLfailpercentcycleT1(:,6),'Color',[1   1  0 ] ,'linewidth',1.5);hold on;
plot(xx,cascLfailpercentcycleT1(:,7),'Color',[1  0   1] ,'linewidth',1.5);hold on;
plot(xx,cascLfailpercentcycleT1(:,8),'Color',[0  1   1] ,'linewidth',1.5);hold on;
plot(xx,cascLfailpercentcycleT1(:,9),'Color',[0.5 0.5 0.5],'linewidth',1.5 );hold on;
plot(xx,cascLfailpercentcycleT1(:,10),'Color',[1 0.5  0] ,'linewidth',1.5);hold on;
plot(xx,cascLfailpercentcycleT1(:,11),'Color',[ 0.67  0   1] ,'linewidth',1.5);
% legend(['absorption index $\beta$=',num2str(Beta(2))],['absorption index $\beta$=',num2str(Beta(3))] ,['absorption index $\beta$=',num2str(Beta(4))],['absorption index $\beta$=',num2str(Beta(5))] ,['absorption index $\beta$=',num2str(Beta(6))],['absorption index $\beta$=',num2str(Beta(7))] ,['absorption index $\beta$=',num2str(Beta(8))],['absorption index $\beta$=',num2str(Beta(9))] ,['absorption index $\beta$=',num2str(Beta(10))] ,['absorption index $\beta$=',num2str(Beta(11))],'Interpreter','latex') ,
legend(['$\beta$=',num2str(Beta(2))],['$\beta$=',num2str(Beta(3))] ,['$\beta$=',num2str(Beta(4))],['$\beta$=',num2str(Beta(5))] ,['$\beta$=',num2str(Beta(6))],['$\beta$=',num2str(Beta(7))] ,['$\beta$=',num2str(Beta(8))],['$\beta$=',num2str(Beta(9))] ,['$\beta$=',num2str(Beta(10))] ,['$\beta$=',num2str(Beta(11))],'Interpreter','latex') ,
xlabel('Key threshold values($\hat\alpha$)','Interpreter','latex');
ylabel({'The complete breakdown ratio ${P_r}$ '},'Interpreter','latex');
set(gca,'FontSize',18); 
% plot(xx,cascLfailpercentcycleT1(:,n),Color(n));
% hols on;
% xlabel('Threshold values($\alpha$)');
% ylabel('The avalanche size P_{f}');


