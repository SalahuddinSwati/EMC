clc;
clear all;


load('MOA_RBF.mat'); % load data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%stream%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stream_data=MOA_RBF; % Data in .mat format, values normalized, last col labels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Important variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global novel;
novel={};
global b_no;
global totaltime;
totaltime=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Important var%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global total_novel_ins;
total_novel_ins=0;
global false_neg_novel_ins;
false_neg_novel_ins=0;
global total_exist_ins;
total_exist_ins=0;
global false_pos_exist_ins;
false_pos_exist_ins=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global num_cluster; %number of clusters
num_cluster=50;
global block_size;
block_size=1000;
init_num_block_train=3; %initial training blocks
global lamda; lamda=0.0009;%%Decat rate
knn=3;
global lofT; lofT=1.5; 
global nooftime; nooftime=3;
global lofk; lofk=10;%%%%%%%%%%%%%%%%%important variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%initial training and test data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    train_data=stream_data(1:block_size*init_num_block_train,:);
    train_data_labels=train_data(:,end);
    
    test_data=stream_data(block_size*init_num_block_train+1:end,1:end-1);
    test_data_labels=stream_data(block_size*init_num_block_train+1:end,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%initial_model_construction%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[train_cls_lb, ia2, cid2] = unique(train_data_labels);

[Model]=initial_model_construction(train_data,num_cluster,train_cls_lb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

correct=0;
fe=0;
test_size=size(test_data,1);
sV=1;
eV=block_size;
j=1;
b_no=1;

ShortMem={};
total_ins=0;
ce=0;
total_ins2=0;
misclall=0;
while eV+block_size<test_size
    %disp(b_no);
    block_data=test_data(sV:eV,:);
    block_lb=test_data_labels(sV:eV,:);
     p = randperm(block_size);
    block_data=block_data(p,:);
    block_lb=block_lb(p);
       
    for i=1: block_size
        total_ins=total_ins+1;
        ex.data=block_data(i,:);
        ex.time=j;
        ex.lb=block_lb(i) ;
        CurrentTime=j;
        clu_centers=cell2mat(Model(:,7));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%all instance classification%%%%%%%%%%%
       
            [ p_label,idx ]=classify(block_data(i,:), Model, knn);
            
            if train_cls_lb(p_label)==block_lb(i)
                
            else
               
                misclall=misclall+1;
            end
        
        %%%%%%%%%%%%%%%%%%%%Existing%%%%%%%%%%%%%%%%%%
        if ismember(block_lb(i),train_cls_lb)
            total_ins2=total_ins2+1;
            [ p_label,idx ]=classify(block_data(i,:), Model, knn);
            
            if train_cls_lb(p_label)==block_lb(i)
                
            else
                ce=ce+1;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [idx, D]=knnsearch(clu_centers,ex.data,'NSMethod','exhaustive');
        r=Model{idx,9};
        if D<= r
            
        else
            
            ShortMem{size(ShortMem,1)+1,1}=ex.data;
            
            ShortMem{size(ShortMem,1),2}=block_lb(i);
            
            ShortMem{size(ShortMem,1),3}=i;
            
        end
        % ex.time
        j=j+1;
    end
    
    
    if (size(ShortMem,1)~=0)

        [ Model, ShortMem, det_new]=detection_novel_class(Model,ShortMem,train_cls_lb, block_lb);
    end
    exist_ins_block=ismember(block_lb,train_cls_lb);
    exist_ins_block= exist_ins_block & ~det_new;
    temp_block_for_classification=block_data(exist_ins_block,:);
    temp_block_label_for_classification=block_lb(exist_ins_block,:);
    for i=1: size(temp_block_for_classification,1)
        [ p_label,idx ]=classify(temp_block_for_classification(i,:), Model, knn);
        
        if train_cls_lb(p_label)==temp_block_label_for_classification(i)
            correct=correct+1;
        else
            fe=fe+1;
            
        end
    end
    [ Model, train_cls_lb ] = update_Model( Model, block_data,block_lb, train_cls_lb,knn,CurrentTime); % Model Update
    sV=eV+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%results%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('block NO = %d \n',b_no);
    cl_Errc=(false_neg_novel_ins+false_pos_exist_ins+fe)*100/total_ins;
    cl_Err=ce*100/total_ins2;
    %%%%%%%%%%%%%%%%%%%%%Classification Error without Novel class%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cl_Err_all=misclall*100/total_ins;
    fprintf('Classification Error without Novel class detection= %f \n',cl_Err_all);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   % fprintf('Classification Error with Novel class detection= %f \n',cl_Err);
   %%%%%%%%%%%%%%%%%%%%%Results with Novel class detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   fprintf('Results with Novel class detection \n');
    fprintf('Combined Error = %f \n',cl_Errc);
    Fn=false_neg_novel_ins*100/total_novel_ins;
    fprintf('novel detection error = %f \n',Fn);
    Fp=false_pos_exist_ins*100/total_exist_ins;
    fprintf('exist detection error = %f \n',Fp);
    fprintf('average time to process block = %f \n',(totaltime/b_no)/60);
     eV=eV+block_size;
     b_no=b_no+1;
end
