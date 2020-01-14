function [Model]=initial_model_construction(data,num_cluster,train_class_labels)

train_data=data(:,1:end-1);
train_data_labels=data(:,end);
%[train_class_labels, ia, cid] = unique(train_data_labels,'stable');
num_replicates=1;
[membership, ctrs] = kmeans(train_data,num_cluster,'Replicates',num_replicates,'Distance','sqEuclidean','MaxIter',1000);
Model={};
for j=1:num_cluster
    clu_pt=find(membership==j);
    
    clu_data=train_data(clu_pt,:);
    clu_labels=train_data_labels(clu_pt,1);
    N_pt=size(clu_data,1);
    clu_center=ctrs(j,:);
    
    LS=sum(clu_data,1);
    SS=sum(clu_data.^2,1);
        clu_r=sqrt(sum(SS/N_pt)-sum((LS/N_pt).^2));
    micro_clu{1,1}=LS; % linear sum
    micro_clu{1,2}=SS; % squared sum
    N_lb=N_pt;%sum(train_data_labels(clu_pt,2));
    N_ub=N_pt-N_lb;
    micro_clu{1,3}=[N_lb N_ub]; % number of labeled and unlabeled instances
    
    LD=[];
    LDC=[];
    RRR=[];
    DCL=[];
    tr_samples={};%for extra point
    perpt=30;%for extra point
    for i=1:size(train_class_labels)
       
        aa=train_data_labels(clu_pt,1)==train_class_labels(i);
         DCL(i)=train_class_labels(i);       
        LD(i)=sum(aa);
      
        rdp=find(aa==1);%%for instances generation

        if LD(i)<=1
           
            tr_samples{i,1}=0;%for instances generation
            tr_samples{i,2}=0;
        else
            X_data=train_data(clu_pt(rdp),:);
            mu=mean(X_data);
            SIGMA=cov(X_data);
        tr_samples{i,1}=mu;%for instances generation
        tr_samples{i,2}=SIGMA;
        end
        %for instances generation
        LDC(i,:)=sum(train_data(clu_pt(aa),:),1);
        cls_LS(i,:)=sum(train_data(clu_pt(aa),:),1);
        cls_SS(i,:)=sum(train_data(clu_pt(aa),:).^2,1);
        
        if LD(i)==0
            RRR(i)=0;
            LDC(i,:)=0;
        else
            RRR(i)=sqrt(sum(cls_SS(i,:)/LD(i))-sum((cls_LS(i,:)/LD(i)).^2));
           
            LDC(i,:)=sum(train_data(clu_pt(aa),:),1)/LD(i);
        end
        
    end
    
    micro_clu{1,4}=LD;
    
    micro_clu{1,5}=0; % time
    
    micro_clu{1,6}=1; % importance
    micro_clu{1,7}=clu_center; 
    micro_clu{1,8}=LDC; % every class center
    micro_clu{1,9}=clu_r; % radius
    micro_clu{1,10}=cls_LS; % class linear sum
    micro_clu{1,11}=cls_SS; % class square sum
    micro_clu{1,12}=RRR; % %every class points radius
    micro_clu{1,13}=1;
    micro_clu{1,14}=DCL;
     micro_clu{1,15}=tr_samples; %%for instances generation
    Model(j,:)=micro_clu;
   
end

 Model  = cal_reability( Model );
end