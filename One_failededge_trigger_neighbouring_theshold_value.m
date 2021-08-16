clc;
clear;
attack=1;
destn=1;
n1=1;
repair=0;
refine=0;
Lfindex=0;
cycle=1; %0==求每个边的theta,1==级联故障。
attackfirst=1;   %循环修复
attackcycle=0; %攻击找边界最大值
t=0.02;
reeverytheta=[];
% B=textread('F:\py\pyscrptis\nx.betweenness_centrality.txt','' , 'headerlines', 0);% 中介度
B=textread('F:\py\Mtb\6_BA1000_3\nx.degree_centrality.txt','' , 'headerlines', 0);% 中介度
% B=textread('F:\py\pyscrptis\nx.degree_centrality.txt','' , 'headerlines', 0);% 中介度
% B=textread('F:\py\pyscrptis\nx.closeness_centrality.txt','' , 'headerlines', 0);% 中介度
B(find(B~=0))=1e+5*(B(find(B~=0))).^(1);
% D=textread('F:\py\pyscrptis\node_degree.txt','', 'headerlines', 0);% DU
D=textread('F:\py\Mtb\6_BA1000_3\node_degree.txt','', 'headerlines', 0);% DU
v=size(B,1);
A=zeros(v,v);
Load=zeros(v,v);
Z1=ones(v,v);
% L=textread('F:\py\pyscrptis\edgelist.txt','', 'headerlines', 0);%%连接关系
L=textread('F:\py\Mtb\6_BA1000_3\edgelist.txt','', 'headerlines', 0);%%连接关系
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
if (attackfirst)   %求每个边的theta
    everytheta=[];
    R=[];
    theta=[];
    attacknumber=0;
    attackeverytheta={};
    recyclerefinenumber=0;
    maxaffectedge=[];
%     retheta=[];
    reeverythetarow=[];
    lrelchangeconnectcell={};
    reeverythetamax=[];
    maxatkline=[];
    for atk=1:1:m
        Loadc=Load;
        attacknumber=attacknumber+1;
        A1=A;
        lchange=[];
        if attack
            i=L(atk,1);
            j=L(atk,2);
%             i=500;
%             j=499;
            %                     random=randi(m);
            %                     i=L(random,1);
            %                     j=L(random,2);
            %             % %攻击开始
            %             i=309;
            %             j=348;
            Linepltofine=[i,j];
            Latk1=[i,j];
            A1(i,j)=0;
            A1(j,i)=0;%i-j边断
            Loadc(i,j)=0;
            Loadc(j,i)=0;
            lchange=[lchange;i;j];
        end
        F0=A1;
        F0((F0~=0))=1;
        
        %是否修复或者改善
        if repair
            i=L0(2);
            j=L0(1);
            F0(i,j)=1;
            F0(j,i)=1;
            lchange=[lchange;i;j];
            %res=multi*max(M(i),M(j));
            res=55;
            A1(i,j)=res;
            A1(j,i)=res;   %i-j
            fprintf('增加一条边%d-%d，增加边的权重为%d', i, j, res)
        end
        if refine
            i=309;
            j=348;
            lchange=[lchange;i;j];
            %i=L(MMAX,1);
            %j=L(MMAX,2);
            F0(i,j)=1;
            F0(j,i)=1;
            %res=multi*[D(i)*D(j)^n1];
            res=1000;
            fprintf('修复一条边%d-%d，修复边的权重为%d\n', i, j, res)
            Z1(i,j)=res/[D(i)*D(j)^n1];
            Z1(j,i)=res/[D(i)*D(j)^n1];
            A1(i,j)=res;
            A1(j,i)=res;   %i-j修复
        end
        A2=A1;
        orlF0=F0;
        F0=orlF0;
        S=sum(A2,2);
        F1=zeros(v,v);
        lchangeconnect=[];
        Lc=sum(Load,2)-sum(Loadc,2);
        nllist=size(lchange,1);
        for j=1:1:nllist
            lconnect=find(L(:,1:2)==lchange(j));
            a=lchange(j);
            nlink=size(lconnect,1);
            for i=1:1:nlink
                if lconnect(i)<=m
                    b=L(lconnect(i),2);
                else
                    b=L(lconnect(i)-m,1);
                end
                isab=ismember([a,b],lchangeconnect);
                            isba=ismember([b,a],lchangeconnect);
                            if ~any(isab(:))||any(isba(:))
                                lchangeconnect=[lchangeconnect;a,b;];
                    if [S(a)*S(b)]~=0
                        F1(a,b)=F0(a,b)*[Z1(a,b)*F2(a,b)*[(B(a)-Lfindex*Lc(a))/S(a)+(B(b)-Lfindex*Lc(b))/S(b)]-1];
                    end
                end
            end
        end
        [F1X,F1Y]=find(F1==max(max(F1)));
        if size(F1X,1)~=1
        maxaffectedge=[maxaffectedge,attacknumber;];    
        end
        theta=[theta;lchange',max(max(F1)),A(lchange(1),lchange(2)),F1X(1),F1Y(1);];
        %triangleF1=triu(F1);
        if max(max(F1))>0
            everytheta=[everytheta,F1(F1~=0)'];
            attackeverytheta{attacknumber}=F1(F1~=0)';
        end
        %&&numel(find(everytheta==F1(a,b)))==0
%         avgtheta=mean(theta,2);
%         maxeverytheta=max(max(everytheta));
%         mineverytheta=min(min(everytheta));
%         avgeverytheta=mean(everytheta,2);
        maxatkline=[sortrows(theta,-3)];
        maxatklinesize=size(theta(theta(:,3)>t),1);
%         maxatkline3length=length(unique(maxatkline(:,3)));
%         maxatklineTsum=sum(maxatkline(:,3)>t);
        maxatklinevar=var(maxatkline(:,3));
        maxatklinestd=std(maxatkline(:,3));
    end
% save('8_BA500_4t01nL1B1degreepredictionedge.mat','maxatkline','maxatklinevar','maxatklinestd');
save('6_BA1000_3t002nL1B1degreepredictionedge.mat','maxatkline','maxatklinevar','maxatklinestd');

end
% save('cascattack8_BA1000_4t01nL02B1degree.mat','cascLfail','cascLfailpercent','LFaillinecell','Lfailcell','erfacell','cascavgerfa','avgcascLfail');
fprintf('结束')