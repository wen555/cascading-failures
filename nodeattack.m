clc;
clear;
attack=1;
destn=1;
n1=0.2;
repair=0;
refine=0;
Lfindex=0;
cycle=1; %0==求每个边的theta,1==级联故障。
attackandrefine=1;
recyclerefine=0;   %循环修复
cascatkrefine=0;   %级联修复
cascattack=0;  %级联故障

attackcycle=0; %攻击找边界最大值
t=0.05;
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
    lconnect=[];
    for atk=1:1:v
        Loadc=Load;
        attacknumber=attacknumber+1;
        A1=A;
        lchange=[];
        llconnect=[];
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
            %Node节点攻击
            [llconnectx,llconnecty]=find(A(atk,:)~=0);
            llconnect= [atk*ones(sum(llconnectx),1),llconnecty'];
            A1((atk-1)*v+llconnecty)=0;
            A1((llconnecty-1)*v+atk)=0;
            Loadc((atk-1)*v+llconnecty)=0;
            Loadc((llconnecty-1)*v+atk)=0;
            lchange=[lchange;llconnect];
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
        Lc=sum(Load,2)-sum(Loadc,2);
        nllist=size(unique(lchange),1);
        lchangeunique=unique(lchange);
        for j=1:1:nllist
       lconnect=[];
           [lconnectx,lconnecty]=find(A(lchangeunique(j),:)~=0);
            a=lchangeunique(j);
            nlink=size(lconnecty,2);
            for i=1:1:nlink
                b=lconnecty(i);
                lconnect=[lconnect;a,b];
                    if [S(a)*S(b)]~=0
                        F1(a,b)=F0(a,b)*[Z1(a,b)*F2(a,b)*[(B(a)-Lfindex*Lc(a))/S(a)+(B(b)-Lfindex*Lc(b))/S(b)]-1];
                    end
            end
        end
%       theta=[theta;lchange',max(max(F1)),A(lchange(1),lchange(2));];
        theta=[theta;max(max(F1)),atk];
        %triangleF1=triu(F1);
        if max(max(F1))>0
            everytheta=[everytheta,F1(F1~=0)'];
            attackeverytheta{attacknumber}=F1(F1~=0)';
        end
        avgtheta=mean(theta,2);
        maxeverytheta=max(max(everytheta));
        mineverytheta=min(min(everytheta));
        avgeverytheta=mean(everytheta,2);
       
%         maxatkline=[sortrows(theta,-3)];
%         maxatklinesize=size(theta(theta(:,3)>t),1);
%         maxatkline3length=length(unique(maxatkline(:,3)));
%         maxatklineTsum=sum(maxatkline(:,3)>t);
%         maxatklinevar=var(maxatkline(:,3));
%         maxatklinestd=std(maxatkline(:,3));
    end
%     save('8_BA1000_4t01nL02B1degree.mat','maxatkline','maxatkline3length','maxatklineTsum','maxatklinevar','maxatklinestd','reeverytheta');
%  save('4_BA100_0t005nL02B1degree.mat','reeverythetarow')
%     [sorttheta atkplace]=sort(theta,'descend');
    %     for line=1:1:size(theta,2);
    %         maxatkline(line,1:4)=[L(atkplace(line),1:2),sorttheta(line),A(L(atkplace(line),1),L(atkplace(line),2))];
    %     end
     [R,P] = corrcoef(B,theta(:,1))
end


if attackcycle
    A3=A;
    erfa=[];
    if destn==1;
        i=L(m,1);
        j=L(m,2);
    elseif destn==2;
        i=feeatk(1,1);
        j=feeatk(1,2);
    else
        random=randi(m);
        i=L(random,1);
        j=L(random,2);
    end
    Latk=[i,j];
    A3(i,j)=0;
    A3(j,i)=0;%i-j边断
    Lfail=1;
    while cycle&&any(any(A3))
        S=sum(A3,2);
        n=0;
        F0=A3;
        F1=zeros(v,v);
        F0((F0~=0))=1;
        for j=1:1:v
            for i=1:1:v
                if (S(j)*S(i))~=0
                    F1(j,i)=F0(j,i)*[Z1(j,i)*F2(j,i)*[B(j)/S(j)+B(i)/S(i)]-1];
                    if (F1(j,i)>t);
                        A3(j,i)=0;
                        n=n+1;
                    end
                end
            end
        end
        erfa=[erfa,max(max(F1))];
        Lfail=[Lfail,n/2];
        if max(max(F1))<=t;
            break;
        end
    end
    avgerfa=mean(erfa,2);
    quantfpercent=sum(Lfail,2)/m;
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
            [Lfaili,Lfailj]=find(triu(F1)>0);
            Lfaillinecell{cyclenumber}=sortrows([Lfaili,Lfailj,F1(triu(F1)>0),A(sub2ind(size(A),Lfaili,Lfailj))],-3);
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
    avgcascLfail=mean(cascLfail,2);
    cascLfailpercent=size(find(cascLfail./m==1),2)/cascnumber;
end
% save('cascattack8_BA1000_4t01nL02B1degree.mat','cascLfail','cascLfailpercent','LFaillinecell','Lfailcell','erfacell','cascavgerfa','avgcascLfail');
fprintf('结束')
