clc;
clear;
cascLfailpercentcycleT=[];
for t=0:0.02:3; 
attack=1;
destn=1;
n1=1;
repair=0;
refine=0;
Lfindex=0;
cycle=1; %0==求每个边的theta,1==级联故障。:0.02:2
   %循环修复
cascatkrefine=0;   %级联修复

cascattack=1;  %级联故障
if cascatkrefine
%     y=load('4_BA100_0t005nL02B1degree.mat','reeverythetarow');
    y=load('8_BA50_4t005nL1B1degreeweak.mat','reeverythetarow');
    reeverythetarow=y.reeverythetarow;
end
attackcycle=0; %攻击找边界最大值
reeverytheta=[];
% B=textread('F:\py\pyscrptis\nx.betweenness_centrality.txt','' , 'headerlines', 0);% 中介度
B=textread('F:\py\Mtb\8_BA50_4\nx.degree_centrality.txt','' , 'headerlines', 0);% 中介度
% B=textread('F:\py\pyscrptis\nx.degree_centrality.txt','' , 'headerlines', 0);% 中介度
% B=textread('F:\py\pyscrptis\nx.closeness_centrality.txt','' , 'headerlines', 0);% 中介度
B(find(B~=0))=1e+5*(B(find(B~=0))).^(1);
% D=textread('F:\py\pyscrptis\node_degree.txt','', 'headerlines', 0);% DU
D=textread('F:\py\Mtb\8_BA50_4\node_degree.txt','', 'headerlines', 0);% DU
v=size(B,1);
A=zeros(v,v);
Load=zeros(v,v);
Z1=ones(v,v);
% L=textread('F:\py\pyscrptis\edgelist.txt','', 'headerlines', 0);%%连接关系
L=textread('F:\py\Mtb\8_BA50_4\edgelist.txt','', 'headerlines', 0);%%连接关系
% 寻找连接边是否存在
[m,n]=size(L);
for i=1:1:m
    A(L(i,1),L(i,2))=[D(L(i,1))*D(L(i,2))]^n1;
end
A=A+A';
%M=mean(A,2);%每一行平均值
S=sum(A,2);
[explorevalue exploresort]=sort(B./S);
BS=B./S;
for i=1:1:m
    BSline(i,1)=BS(L(i,1));
    BSline(i,2)=BS(L(i,2));
end

exBSorder=[explorevalue exploresort];
exBSReverseorder=flipud(exBSorder);
for i=1:1:v
    for j=1:1:v
        if S(j)~=0
            exAS(i,j)=A(i,j)/S(j);
        end
    end
end
[exASi,exASj]=find(exAS~=0);
exASorder=sortrows([exASi,exASj,exAS(exAS~=0),A(sub2ind(size(A),exASi,exASj))],1);

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
%         cascatki=3;
%         cascatkj=15;

        Lcascatk=[cascatki,cascatkj];
        A2(cascatki,cascatkj)=0;
        A2(cascatkj,cascatki)=0;
        if cascatkrefine&&ismember(Lcascatk,reeverythetarow(:,1:2),'rows')
            cascam=cascam+1;
            [cascatkjnumber,cascatkjrows]=ismember(Lcascatk,reeverythetarow(:,1:2),'rows');
            %half1
            %生成一个仅包含[0,1]的数组



%             p1 = unifrnd(0,1,1,m);%(0,1)均匀分布中随机抽取一些数
%                 if p1(1,cascatk)>=0.5
%                     p1(1,cascatk) = 1;
%                 else p1(1,cascatk) = 0;
%                 end;
%         p1=[];
        p1(1,cascatk)=0;
     if p1(1,cascatk)==1
            cascatkires=reeverythetarow(cascatkjrows,3);
            cascatkjres=reeverythetarow(cascatkjrows,4);
            res=(D(cascatkires)*D(cascatkjres))^n1;
            A2(cascatkires,cascatkjres)=res;
            A2(cascatkjres,cascatkires)=res;
%           修复两边
            %half2
     elseif p1(1,cascatk)==2
            cascatkires=reeverythetarow(cascatkjrows+1,3);
            cascatkjres=reeverythetarow(cascatkjrows+1,4);
            res=(D(cascatkires)*D(cascatkjres))^n1;
            A2(cascatkires,cascatkjres)=res;
            A2(cascatkjres,cascatkires)=res;
     else p1(1,cascatk)==3
           cascatkires=reeverythetarow(cascatkjrows,3);
            cascatkjres=reeverythetarow(cascatkjrows,4);
            res=(D(cascatkires)*D(cascatkjres))^n1;
            A2(cascatkires,cascatkjres)=res;
            A2(cascatkjres,cascatkires)=res;
            cascatkires=reeverythetarow(cascatkjrows+1,3);
            cascatkjres=reeverythetarow(cascatkjrows+1,4);
            res=(D(cascatkires)*D(cascatkjres))^n1;
            A2(cascatkires,cascatkjres)=res;
            A2(cascatkjres,cascatkires)=res;
     end
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
            [Lfaili,Lfailj]=find(triu(F1)~=0);
            Lfaillinecell{cyclenumber}=sortrows([Lfaili,Lfailj,F1(triu(F1)~=0),A(sub2ind(size(A),Lfaili,Lfailj))],-3);
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
    avgcascLfail=(mean(cascLfail,2)-1)/m;
    cascLfailpercent=size(find(cascLfail./m>=1),2)/cascnumber;
end
cascLfailpercentcycleT=[cascLfailpercentcycleT;t,avgcascLfail];
end
% cascLfailpercentcycleT=cascLfailpercentcycleT';&~(ismember([i,j],[cascatkires,cascatkjres],'rows'))&~(ismember([i,j],[cascatkjres,cascatkires],'rows'))
save('8_BA50_4t005nL1B1degreecascadingfailNorefine.mat','cascLfailpercentcycleT')
% save('8_BA50_4t005nL1B1degreecascadingfailThalf1weak.mat','cascLfailpercentcycleT')
% save('8_BA50_4t005nL1B1degreecascadingfailThalf2strong.mat','cascLfailpercentcycleT')
% save('8_BA50_4t005nL1B1degreecascadingReTstrong.mat','cascLfailpercentcycleT')
fprintf('结束')