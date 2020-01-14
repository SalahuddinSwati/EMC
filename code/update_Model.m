function [ Model, tr_cls_lb ] = update_Model( Model, data,block_labels,tr_cls_lb,knn,CurrentTime)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global lamda;
[~, ib] = ismember(block_labels, tr_cls_lb);
temp_data=data(ib~=0,:);
new_cls_data=data(ib==0,:);
temp_lb=block_labels(ib~=0);
new_cls_lbs=block_labels(ib==0);
new_label=[];
if size(new_cls_lbs,1)>1
    new_label=new_cls_lbs(1);
end
for i=1:size(data,1)
    if ismember(block_labels(i),tr_cls_lb)
        [ p_label,idx ]=classify(data(i,:), Model, knn);
        
        if tr_cls_lb(p_label)==block_labels(i)
            
            Model(idx,6)=num2cell(cell2mat(Model(idx,6))+1);
            
        else
            
            Model(idx,6)=num2cell(cell2mat(Model(idx,6))-1);
        end
    end
end

[ Model ] = break_MC( Model,CurrentTime );
impr=cell2mat(Model(:,6));
impr=impr.*2.^(-lamda.*(CurrentTime-cell2mat(Model(:,13))));
Model(:,6)=num2cell(impr);


impr=cell2mat(Model(:,6));

Model(impr<0.1,:)=[];
z=1;
y=1;
rem=[];
false_ins=[];
for i=1:size(temp_data,1)
    ex.data=temp_data(i,:);
    ex.lb=temp_lb(i);
    clu_centers=cell2mat(Model(:,7));
    [idx, D]=knnsearch(clu_centers,ex.data,'NSMethod','exhaustive');
    r=Model{idx,9};
    if D<=r
        
        [ Model ] = update_micro(Model,ex,idx, tr_cls_lb, CurrentTime);
        
        rem(z)=i;
        z=z+1;
    end
end

if size(rem,1)>0
    temp_data(rem,:)=[];
    temp_lb(rem)=[];
end
if size(new_cls_data,1)>1
    temp_data=[temp_data;new_cls_data];
    temp_lb=[temp_lb; new_cls_lbs];
end
global num_cluster;
global block_size;
num_replicates=5;

T=sum(cell2mat(Model(:,4)),1);
idx=find(T==0);
if size(idx,2)~=0
    for i=1:size(Model,1)
        Model{i,4}(idx)=[];
        Model{i,8}(idx,:)=[];
        Model{i,10}(idx,:)=[];
        Model{i,11}(idx,:)=[];
        Model{i,12}(idx)=[];
        Model{i,14}(idx)=[];
        Model{i,15}(idx,:)=[];
    end
    tr_cls_lb(idx)=[];
end
if size(new_label,1)>=1
    
    for i=1:size(Model,1)
        Model{i,4}(end+1)=0;
        
        Model{i,8}(end+1,:)=0;
        Model{i,10}(end+1,:)=0;
        Model{i,11}(end+1,:)=0;
        Model{i,12}(end+1)=0;
        Model{i,14}(end+1)=new_label;
        
        Model{i,15}{end+1,1}=0;
        Model{i,15}{end,2}=0;
        
    end
    tr_cls_lb=[tr_cls_lb;new_label];
end
K=ceil(num_cluster*size(temp_data,1)/block_size);
if K>=1
    [membership, ctrs] = kmeans(temp_data,K,'Replicates',num_replicates,'Distance','sqEuclidean','MaxIter',1000);
    tr_samples={};%for extra point
    perpt=30;%for extra point
    
    for j=1:K
        
        clu_pt=find(membership==j);
        
        clu_data=temp_data(clu_pt,:);
        N_pt=size(clu_data,1);
%         if N_pt==1
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             ex.data=clu_data;
%             ex.lb=temp_lb(clu_pt);
%             clu_centers=cell2mat(Model(:,7));
%             [idx, D]=knnsearch(clu_centers,ex.data,'NSMethod','exhaustive');
%             
%             [ Model ] = update_micro(Model,ex,idx, tr_cls_lb, CurrentTime);
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%
%         else
            clu_center=ctrs(j,:);
            
            LS=sum(clu_data,1);
            SS=sum(clu_data.^2,1);
            
            clu_r=sqrt(sum(SS/N_pt)-sum((LS/N_pt).^2));
            micro_clu{1,1}=LS;
            micro_clu{1,2}=SS;
            
            N_lb=N_pt;%sum(train_data_labels(clu_pt));
            N_ub=N_pt-N_lb;
            micro_clu{1,3}=[N_lb N_ub];
            cls_LS=[];
            cls_SS=[];
            
            RRR=[];
            LD=[];
            LDC=[];
            DCL=[];
            for i=1:size(tr_cls_lb)
                aa=temp_lb(clu_pt)==tr_cls_lb(i);
                DCL(i)=tr_cls_lb(i);
                LD(i)=sum(aa);
                
                %for extra point
                rdp=find(aa==1);%for extra point
                sz_rp=size(rdp,1);%for extra point
                red=randperm(sz_rp);%for extra point
                if LD(i)<=1
                    
                    tr_samples{i,1}=0;%for extra point
                    tr_samples{i,2}=0;
                else
                    X_data=temp_data(clu_pt(rdp),:);
                    mu=mean(X_data);
                    SIGMA=cov(X_data);
                    tr_samples{i,1}=mu;%for extra point
                    tr_samples{i,2}=SIGMA;
                end
                %for extra point
                LDC(i,:)=sum(temp_data(clu_pt(aa),:),1)/LD(i);
                cls_LS(i,:)=sum(temp_data(clu_pt(aa),:),1);
                cls_SS(i,:)=sum(temp_data(clu_pt(aa),:).^2,1);
                if LD(i)==0
                    RRR(i)=0;
                    LDC(i,:)=0;
                else
                    RRR(i)=sqrt(sum(cls_SS(i,:)/LD(i))-sum((cls_LS(i,:)/LD(i)).^2));
                    LDC(i,:)=sum(temp_data(clu_pt(aa),:),1)/LD(i);
                end
                
            end
            micro_clu{1,4}=LD;
            
            micro_clu{1,5}=0;
            
            micro_clu{1,6}=1;
            micro_clu{1,7}=clu_center;
            micro_clu{1,8}=LDC;
            micro_clu{1,9}=clu_r;
            micro_clu{1,10}=cls_LS;
            micro_clu{1,11}=cls_SS;
            micro_clu{1,12}=RRR;
            micro_clu{1,13}=CurrentTime;
            micro_clu{1,14}=DCL;
            micro_clu{1,15}=tr_samples; % for extra points
            m_idx=size(Model,1)+1;
            Model(m_idx,:)=micro_clu;
            
            
%         end
    end
end

[ Model ] = cal_reability( Model );
end
%     im=1;
%     R=[];
%     for i=1:1000:25000
%         R=[R im*2.^(-0.00003*i)];
%         im=im*2.^(-0.00003*i);
%     end
%     plot(R)