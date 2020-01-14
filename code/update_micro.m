function [ Model ] = update_micro(Model,ex,idx, tr_cls, CurrentTime)
   
    lbd=Model{idx,3};
    
        lbd(1,1)=lbd(1,1)+1;
    
    LS=Model{idx,1}+ex.data;
    SS=Model{idx,2}+ex.data.^2;
    N_pt=sum(Model{idx,3})+1;
    clu_center=LS/N_pt;
    clu_r=sqrt(sum(SS/N_pt)-sum((LS/N_pt).^2));

    Model{idx,1}=LS;
    Model{idx,2}=SS;
    Model{idx,3}=lbd;
    
        tr_cls_pt=find(tr_cls==ex.lb);
        Model{idx,4}(1,tr_cls_pt)=Model{idx,4}(1,tr_cls_pt)+1;
        LSC=Model{idx,10}(tr_cls_pt,:)+ex.data;
        Model{idx,10}(tr_cls_pt,:)=LSC;
        SSC=Model{idx,11}(tr_cls_pt,:)+ex.data.^2;
        Model{idx,11}(tr_cls_pt,:)=SSC;
        nptc=Model{idx,4}(1,tr_cls_pt);
        Model{idx,8}(tr_cls_pt,:)=LSC/nptc;
        Model{idx,12}(1,tr_cls_pt)=sqrt(sum(SSC/nptc)-sum((LSC/nptc).^2));
        Model{idx,13}=CurrentTime;
    Model{idx,7}=clu_center;
    Model{idx,9}=clu_r;
end

