    function [ Model, ShortMem, det_new] = detection_novel_class(Model,ShortMem,tr_cls, block_lb)
    global num_cluster;
    global block_size;
    K=ceil(num_cluster*size(ShortMem,1)/block_size);
    num_replicates=1;
    class_center=cell2mat(Model(:,8));
    M_Clu_C=cell2mat(Model(:,7));
    ShortMem_data=cell2mat(ShortMem(:,1));
    cls_dim=size(ShortMem_data,2);
    sh_lb=cell2mat(ShortMem(:,2));
    sh_lb=sh_lb(:,1);
    Total_new=~ismember(block_lb,tr_cls);
    DimSz=size(Model{1,1},2);


    global lofT; global nooftime; global lofk; %%%%%%%%%%%%%%%%%important variable


    global novel;

    global b_no;

    global total_novel_ins;

    global false_neg_novel_ins;

    global total_exist_ins;

    global false_pos_exist_ins;

    %     total_new_pt=[];

    [membership, ctrs] = kmeans(ShortMem_data,K,'Replicates',num_replicates,'Distance','sqEuclidean');
    det_new=zeros(size(block_lb,1),1);
    global totaltime;
    lof_data={};
    cov_lof_data={};
    tic
    for i=1: size(tr_cls,1)
        lof_data_temp=[];
        cls_data = cellfun(@(x) x(i,:),Model(:,15),'UniformOutput',false);
        cls_num = cellfun(@(y) y(i),Model(:,4),'UniformOutput',false);
        cls_cntr = cellfun(@(z) z(i,:),Model(:,8),'UniformOutput',false);
        cls_cntr=cell2mat(cls_cntr);
        cls_num=cell2mat(cls_num);
        indx=find(cls_num~=0);
        cls_num=cls_num(indx);

        cls_data=cls_data(indx);
        cls_cntr=cls_cntr(indx,:);

        clsDsz=size(cls_data,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        totalsum=sum(cls_num);

        clsDim=size(cls_cntr,2); %change
        totalex=nooftime*clsDim;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for m=1:clsDsz
            X_data=[];
            if cls_num(m)<=2
                X_data= cls_cntr(m,:);
            else
                %%%%%%%%%%%%%%%%%
                numofsample=ceil((totalex*cls_num(m))/totalsum);
                %%%%%%%%%%%%%%
                if size(cls_data{m}{1,1},2)>1
                    mu=cls_data{m}{1,1};
                    sigma=cls_data{m}{1,2};

                    X_data= mvnrnd(mu,sigma,numofsample);
                end
            end
            %         if size(X_data,2)
            %             disp(X_data);
            %         else
            lof_data_temp=[lof_data_temp' X_data']';
            %         end


        end

        lof_data{i}=lof_data_temp;
       
        if size(lof_data_temp,1)>cls_dim
             cov_lof_data{i}=nearestSPD(cov(lof_data_temp));
            [k_index, k_dist] = knnsearch(lof_data_temp,lof_data_temp,'k',lofk+1,'nsmethod','exhaustive','Distance','mahalanobis', 'IncludeTies',true, 'Cov',nearestSPD(cov(lof_data_temp)));
            %Ignore first element(itself) at nearest neighbors
            k_index = cellfun(@(x) x(2:end),k_index,'UniformOutput',false);
            numneigh = cellfun('length',k_index);
            %Get k-distance
            k_dist1 = cell2mat(cellfun(@(x) x(end),k_dist,'UniformOutput',false));
            k_indexcls{i}=k_index;
            numneighcls{i}=numneigh;
            k_dist1cls{i}=k_dist1;
        else
             cov_lof_data{i}=0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%

    time1=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    for j=1:size(ShortMem_data,1)
        strt=tic;
        xx=ShortMem_data(j,:);


        cls_lof=[];

        for i=1: size(tr_cls,1)
            f=0;
            lof_dt=lof_data{i};
            cov_dt=cov_lof_data{i};
            if size(lof_dt,1)>clsDim


                [k_indexx, k_distx] = knnsearch(lof_dt,xx,'k',lofk,'nsmethod','exhaustive','Distance','mahalanobis', 'IncludeTies',true, 'Cov',cov_dt);
                %Ignore first element(itself) at nearest neighbors
                %k_index = cellfun(@(x) x(2:end),k_index,'UniformOutput',false);
                numneighx = cellfun('length',k_indexx);
                %Get k-distance
                k_dist1x = cell2mat(cellfun(@(x) x(end),k_distx,'UniformOutput',false));

                T_k_indx= k_indexcls{i};
                T_k_indx{end+1}=cell2mat(k_indexx);

                T_numneigh=numneighcls{i};
                T_numneigh(end+1)=numneighx;

                T_k_dist1=k_dist1cls{i};
                T_k_dist1(end+1)=k_dist1x;

                lof_dt=[lof_dt;xx];


                [suspicious_index lof] = LOF_maha(lof_dt, lofk,T_k_indx,T_numneigh, T_k_dist1,cov_dt);


                cls_lof(i)=lof;
                if lof<lofT
                    f=1;
                    break;
                end
            else

                cls_lof(i)=0;
            end
        end

        n_sz=size(novel,1)+1;

        if f==0
            novel{n_sz,1}=1;
            %
            det_new(ShortMem{j,3})=1;

        else
            %
            novel{n_sz,1}=0;
        end
        novel{n_sz,2}=cls_lof;

        novel{n_sz,3}=b_no;
        tt=toc(strt);
        time1=time1+tt;
        %totaltime=totaltime+time1;
       % fprintf('Time for %d / %d = %f \n',j, size(ShortMem_data,1),tt/60);
    end
    %%
    if sum(det_new)< 50
        det_new=0;
    end
      time2=toc;
    totaltime=totaltime+time2;
      total_exist_ins=total_exist_ins+(block_size-sum(Total_new));
   
    total_novel_ins=total_novel_ins+sum(Total_new);
      false_pos_exist_ins=false_pos_exist_ins+sum(~Total_new & det_new);
     false_neg_novel_ins=false_neg_novel_ins+sum(Total_new & ~det_new);
    ShortMem={};
    end

