function [ p_label, idx] = classify( ex, Model, knn )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

  clu_centers=cell2mat(Model(:,7));
  no_of_cls=length(cell2mat(Model(1,4)));
[idx, D]=knnsearch(clu_centers,ex,'NSMethod','exhaustive','k',knn);
        dist=zeros(no_of_cls,1);
        for l=1:no_of_cls
            for j=1:knn
                ix=idx(j);
                LC=cell2mat(Model(ix,8));
                R=cell2mat(Model(ix,5));
                pr=cell2mat(Model(ix,4));
                dist(l)=dist(l)+(R*pr(l)/sum(pr))/norm(ex-LC(l,:));
               
            end
            
        end
        
       [v, p_label]=max(dist);
end

