clc;
clear;
attack=1;
destn=2;
n1=1;
repair=0;
refine=0;
Lfindex=0;
cycle=1; %0==求每个边的theta,1==级联故障。4_SW500_4_08
attackcycle=1; %攻击找边界最大值8_BA500_4 8_SW500_8_05
t=0.1;
edgenotzero=[];
reeverytheta=[];
% B=textread('F:\py\pyscrptis\nx.betweenness_centrality.txt','' , 'headerlines', 0);% 中介度
B=textread('F:\py\Mtb\8_SW500_8_05\nx.degree_centrality.txt','' , 'headerlines', 0);% 中介度
% B=textread('F:\py\pyscrptis\nx.degree_centrality.txt','' , 'headerlines', 0);% 中介度
% B=textread('F:\py\pyscrptis\nx.closeness_centrality.txt','' , 'headerlines', 0);% 中介度
B(find(B~=0))=1e+5*(B(find(B~=0))).^(1);
% D=textread('F:\py\pyscrptis\node_degree.txt','', 'headerlines', 0);% DU
D=textread('F:\py\Mtb\8_SW500_8_05\node_degree.txt','', 'headerlines', 0);% DU
v=size(B,1);
A=zeros(v,v);
Load=zeros(v,v);
Z1=ones(v,v);
% L=textread('F:\py\pyscrptis\edgelist.txt','', 'headerlines', 0);%%连接关系
L=textread('F:\py\Mtb\8_SW500_8_05\edgelist.txt','', 'headerlines', 0);%%连接关系
% 寻找连接边是否存在
[m,n]=size(L);
for i=1:1:m
    A(L(i,1),L(i,2))=[D(L(i,1))*D(L(i,2))]^n1;
end
A=A+A';
%M=mean(A,2);%每一行平均值
S=sum(A,2);
Ss=S./sum(A~=0,2);
[explorevalue exploresort]=sort(B./S);
BS=B./S;
for i=1:1:m
    BSline(i,1)=BS(L(i,1));
    BSline(i,2)=BS(L(i,2));
end

exBSorder=[explorevalue exploresort];
exBSReverseorder=flipud(exBSorder);
exAS=zeros(v);
% for i=1:1:v
%     for j=1:1:v
%         if S(j)~=0
%             exAS(i,j)=A(i,j)/S(j);
%         end
%     end
% end
for i=1:1:m
        if Ss(L(i,2))~=0
            exAS(L(i,1),L(i,2))=(A(L(i,1),L(i,2))/(Ss(L(i,2))));
        end
end
exAS=exAS+exAS';


[exASi,exASj]=find(exAS~=0);
% exASorder=sortrows([exASi,exASj,exAS(exAS~=0)],1);
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
if attackcycle
    A3=A;
    erfa=[];
    if destn==1
        
%         i=5;
%         j=96;
        i=5;
        j=317;
    elseif destn==2
        random=randi(m);
        i=L(random,1);
        j=L(random,2);
    end
    Latk=[i,j];
    A3(i,j)=0;
    A3(j,i)=0;%i-j边断
    linerefine=[];
    if refine
        i=1;
        j=38;
    linerefine=[i,j];
    A3(i,j)=[D(i)*D(j)]^n1;
    A3(j,i)=[D(i)*D(j)]^n1;
%         i=26;
%         j=21;
%     linerefine=[linerefine;[i,j]];
%     A3(i,j)=[D(i)*D(j)]^n1;
%     A3(j,i)=[D(i)*D(j)]^n1;
    end
    
    Lfail=1;
    cyclenumber=0;
    while cycle&&any(any(A3))
        S=sum(A3,2);
        n=0;
        F0=A3;
        edgenotzero1=[];
        cyclenumber=cyclenumber+1;
        F1=zeros(v,v);
        F0((F0~=0))=1;
        for j=1:1:v
            for i=1:1:v
                if (S(j)*S(i))~=0
                    F1(j,i)=F0(j,i)*[Z1(j,i)*F2(j,i)*[B(j)/S(j)+B(i)/S(i)]-1];
                    if (F1(j,i)>t);
                     A3(j,i)=0;
                        n=n+1;
                        edgenotzero1=[edgenotzero1;Latk,cyclenumber,n,F1(j,i),i,j,Load(i,j),(1+F1(j,i))*Load(i,j),exAS(i,j)] ;  
                    end
                end
            end
        end
        if any(edgenotzero1)
        edgenotzero=[edgenotzero;sortrows(edgenotzero1,-5)] ;
        end
        erfa=[erfa,max(max(F1))];
        Lfail=[Lfail,n/2];
        if max(max(F1))<=t
            break;
        end
    end
    if edgenotzero
    edgenotzero=edgenotzero(2:2:end,:);  
    edgenotzero=sortrows(edgenotzero,6);
avgerfa=mean(erfa,2);
    quantfpercent=sum(Lfail,2)/m;
geterfa=edgenotzero(:,9);
Ci=edgenotzero(:,8);
weightratio=edgenotzero(:,10);
    end
end
txtfile( ['F:\py\Mtb\8_SW500_8_05\','edgenotzero.txt'],edgenotzero)
% save('cascattack8_BA1000_4t01nL02B1degree.mat','cascLfail','cascLfailpercent','LFaillinecell','Lfailcell','erfacell','cascavgerfa','avgcascLfail');
if edgenotzero
figure(1)
subplot(3,1,1)
plot(Ci,geterfa,'.')
subplot(3,1,2)
plot(Ci,weightratio,'.')
subplot(3,1,3)
plot(weightratio,geterfa,'.')
figure(2)
plot(edgenotzero(:,5))
fprintf('结束')


 figure(4);
 subplot(4,2,1)
 plot(edgenotzero(:,10))
xlabel('Failed edge with time sequence ')
ylabel('Edge weight ratio')

 subplot(4,2,2)
 plot(edgenotzero(:,8))
 ylabel('Initial load')
xlabel('Failed edge with time sequence')
 subplot(4,2,3)
  plot(edgenotzero(:,9))
  ylabel('Updated load ')
xlabel('Failed edge with time sequence')
   subplot(4,2,4)
     plot(edgenotzero(:,5));
     ylabel('The updated threshold value')
xlabel('Failed edge with time sequence')
end

     