clc;
clear;
attack=1;
destn=1;
n1=1;
repair=0;
refine=0;
Lfindex=0;
cycle=1; %0==求每个边的theta,1==级联故障。
attackandrefine=1;
recyclerefine=1;   %循环修复
t=0.05;
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
if (attackandrefine)   %求每个边的theta
    everytheta=[];
    R=[];
    theta=[];
    attacknumber=0;
    attackeverytheta={};
    recyclerefinenumber=0;
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
%             i=3;
%             j=15;
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
        if recyclerefine
            rethetaiLafirst=[];
            retheta=[];
            reeverythetaiL=[];
            reFOmax=[];
            reeverythetanumberj=[];
            for numberi=1:1:size(Linepltofine,2);
                iL=Linepltofine(numberi);
                connect=find(L(:,1:2)==iL);
                iconnect=[];
                relchange=[];
                for i1=1:1:size(connect,1)
                    if connect(i1)<=m
                        b=L(connect(i1),2);
                    else
                        b=L(connect(i1)-m,1);
                    end
                    iconnect=[iconnect;iL,b;];
                end
                for jj=1:1:size(iconnect,1)
                    numberj=iconnect(jj,2);
                    if (ismember(numberj,Linepltofine))
                         continue;
                    end
                    % (ismember(numberj,iconnect))
                    j=numberj;
                    if numberi==1;
                    nn=2;
                    else
                    nn=1;
                    end 
                    iL=Linepltofine(nn);
                    if (ismember([j,iL],L(:,1:2),'rows')||(ismember([iL,j],L(:,1:2),'rows')))
                        continue;
                    end
                    recyclerefinenumber=recyclerefinenumber+1;
                    A1=A2;
                    reF0=F0;
                    reF0(iL,j)=1;
                    reF0(j,iL)=1;
                    rethetaiLa=[];
                    relchange=lchange;
                    relchange=[relchange;j];
                    res=(D(iL)*D(j))^n1;
%                   res=1;
                    A1(iL,j)=res;
                    A1(j,iL)=res; 
                    %i-j
                    %fprintf('增加一条边%d-%d，增加边的权重为%d', i, j, res)
                    S=sum(A1,2);
                    reF1=zeros(v,v);
                    lrelchangeconnect=[];
                    nllist=size(relchange,1);
                    for j=1:1:nllist
                        lconnect=find(L(:,1:2)==relchange(j));
                        a=relchange(j);
                        nlink=size(lconnect,1);
                        for i=1:1:nlink
                            if lconnect(i)<=m
                                b=L(lconnect(i),2);
                            else
                                b=L(lconnect(i)-m,1);
                            end
                            if ismember([a,b],Linepltofine,'rows')||ismember([b,a],Linepltofine,'rows')
                                continue;
                            end
                            isab=ismember([a,b],lrelchangeconnect);
                            isba=ismember([b,a],lrelchangeconnect);
                            if ~any(isab(:))||any(isba(:))
                                lrelchangeconnect=[lrelchangeconnect;a,b;];
                                if [S(a)*S(b)]~=0
                                    reF1(a,b)=reF0(a,b)*[Z1(a,b)*F2(a,b)*[B(a)/S(a)+B(b)/S(b)]-1];
                                    if reF1(a,b)~=0
                                        retheta=[retheta;Linepltofine,iL,numberj,a,b,reF1(a,b);];%iL,numberj 修复的边，a,b与之相关的边；
                                    end
                                    if iL==a;
                                       rethetaiLa=[rethetaiLa;Linepltofine,iL,numberj,a,b,reF1(a,b),A(iL,numberj);]  ;
                                    end
                                end
                            end
                        end
                    end
              
                    rethetaiLa=sortrows(rethetaiLa,-7);%从大到小
          
                    rethetaiLafirst=[rethetaiLafirst;rethetaiLa(1,:)];  
                  
%                     rethetaiLafirst=sortrows(rethetaiLafirst,7);%排列从小到大，去修复后第一行 strong 修复
                    rethetaiLafirst=sortrows(rethetaiLafirst,7); %强修复修复
                    lrelchangeconnectcell{recyclerefinenumber}=lrelchangeconnect;
                    if max(max(reF1))>0
                        [xp,yp]=find(reF1==max(max(reF1)));
                        reeverytheta=[reeverytheta;Linepltofine,max(max(reF1)),iL,numberj,A(Linepltofine(1),Linepltofine(2)),A(iL,numberj);];
%                         sortrows(retheta(ffre,:),7)
                        reFOmax=[reFOmax,max(max(reF1))];
                    end
%                     reavgtheta=mean(retheta(:,7),2);
%                     remaxtheta=max(max(retheta(:,7)));
%                     remintheta=min(min(retheta(:,7)));
                    reavgeverytheta=mean(reeverytheta(:,3),2);
                     
                end
            end
            
            reeverytheta=sortrows(reeverytheta,-3);
            reeverythetaiL=[reeverythetaiL;Linepltofine,Linepltofine(1),Linepltofine(1)];
            reeverythetanumberj=[reeverythetanumberj;Linepltofine,Linepltofine(2),Linepltofine(2)];
            [resum,rerows]=ismember(reeverythetaiL,rethetaiLafirst(:,[1 2 3 5]),'rows');
            [resum1,rerows1]=ismember(reeverythetanumberj,rethetaiLafirst(:,[1 2 3 5]),'rows');
        if resum&resum1
            reeverythetarow=[reeverythetarow;rethetaiLafirst(rerows,:);rethetaiLafirst(rerows1,:)];
        end
            
        end
    end
%save('8_BA1000_4t01nL02B1degree.mat','maxatkline','maxatkline3length','maxatklineTsum','maxatklinevar','maxatklinestd','reeverytheta');
 save('8_BA50_4t005nL1B1degreeweak.mat','reeverythetarow')
%[sorttheta atkplace]=sort(theta,'descend');
    %     for line=1:1:size(theta,2);
    %         maxatkline(line,1:4)=[L(atkplace(line),1:2),sorttheta(line),A(L(atkplace(line),1),L(atkplace(line),2))];
    %     end
end
fprintf('结束')