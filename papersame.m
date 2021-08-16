clc;
clear;
cascLfailpercentcycleT=[];
cascLfailpercentfailavgT=[];
for t=0.01
attack=1;
destn=1;
n1=0.1;
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
B=textread('F:\py\Mtb\4_BA100_2\nx.degree_centrality.txt','' , 'headerlines', 0);% 中介度
B(find(B~=0))=1e+5*(B(find(B~=0))).^(1);
% D=textread('F:\py\pyscrptis\node_degree.txt','', 'headerlines', 0);% DU
D=textread('F:\py\Mtb\4_BA100_2\node_degree.txt','', 'headerlines', 0);% DU
v=size(B,1);
A=zeros(v,v);
Load=zeros(v,v);
Z1=ones(v,v);
% L=textread('F:\py\pyscrptis\edgelist.txt','', 'headerlines', 0);%%连接关系
L=textread('F:\py\Mtb\4_BA100_2\edgelist.txt','', 'headerlines', 0);%%连接关系
% 寻找连接边是否存在
[m,n]=size(L);
for i=1:1:m
    A(L(i,1),L(i,2))=[D(L(i,1))*D(L(i,2))]^n1;
end
A=A+A';
%M=mean(A,2);%每一行平均值
S1=sum(A,2);
 A2=A;
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

% for j=1:1:v
%     if S(j)~=0
%         for i=1:1:v
%             Load(j,i)=A(j,i)/S(j)*B(j);
%         end
%     end
% end
% Loadc=Load;
% F2=zeros(v,v);
% for j=1:1:v
%     for i=1:1:v
%         if (B(i)*S(j)+B(j)*S(i))~=0
%             F2(j,i)=S(i)*S(j)/(B(i)*S(j)+B(j)*S(i));
%         end
%     end
% end
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
        Lcas=[];
        A2=A;
        cascatki=L(cascatk,1);
        cascatkj=L(cascatk,2);
        Lcascatk=[cascatki,cascatkj];
        Lcas=Lcascatk;
%         A(cascatki,cascatkj)
        A2(cascatki,cascatkj)=0;
        A2(cascatkj,cascatki)=0;      
        Lfail=1;
        F0=zeros(v,v);
        Lfailline=[];
        cyclenumber=0;
        Lfaillinecell={};
        cascnumber=cascnumber+1; 
        while cycle&any(any(A2))
            S=sum(A2,2);
            A3=A2;
            F0(find(A2~=0))=1;
            n=0;
            F1=zeros(v,v);
            cyclenumber=cyclenumber+1;
            for j=1:1:v
                for i=1:1:v
                    if any(Lcas==i)||any(Lcas==j) 
                    if F0(i,j)*(S(i)+S(j))*A(j,i)~=0
                        F1(j,i)=F0(j,i)*([A(j,i)*(S(i)+S(j))^(-1)+A3(j,i)]/A(j,i)-1);
                      if (F1(j,i)>t);
                            A2(j,i)=0;
                            Lcas=[Lcas,i,j];
                            n=n+1;
                      else
                           A2(j,i)=A(j,i)*(S(i)+S(j))^(-1)+A3(j,i);
                      end
                    end
                    end
                end
                as=size(Lcas,2);
                 if size(Lcas,1)==0;
                break;
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
    avgcascLfail=mean(cascLfail,2)./m;
    cascLfailpercent=size(find(cascLfail./m==1),2)/cascnumber;  
end
cascLfailpercentcycleT=[cascLfailpercentcycleT;t,cascLfailpercent];
cascLfailpercentfailavgT=[cascLfailpercentfailavgT;t,avgcascLfail];
end
% cascLfailpercentcycleT=cascLfailpercentcycleT';
% save('6_BA500_3nL1B1degreecascadingfailT.mat','cascLfailpercentcycleT')
fprintf('over')
    

