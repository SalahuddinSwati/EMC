function [ Model ] = cal_reability( Model )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    lb_data=cell2mat(Model(:,4));
    mx=max(lb_data,[],2);
    sm=sum(lb_data,2);
   purity=max(lb_data,[],2)./sum(lb_data,2);
   Model(:,5)=num2cell(purity);
%     H=0;
%     PD=[];
%     total_lb=sum(sum(lb_data));
%     cls_lb=sum(lb_data,1);
%     for i=1:size(lb_data,2)
% 
%         H=H+(-cls_lb(i)/total_lb)*log(cls_lb(i)/total_lb);
%         PD(i)=cls_lb(i)/total_lb;
%     end
%         
%     for i=1:size(lb_data,1)
%         HC=0;
%         vE=0;
%         if sum(lb_data(i,:))==0
%             Model{i,5}=0;
%         else
%         for j=1:size(lb_data,2)
%             PC=lb_data(i,j)/sum(lb_data(i,:));
%             if PC~=0
%               HC=HC+(-PC)*log(PC);
%                
%             end
%             vE=vE+(PC-PD(j))/PD(j);
%         end
%             
%             CR=(H-HC)/H;
%             
%             CP=1/(1+exp(-vE));
%             R=CR*CP;
%             if isnan(R)
%                 Model{i,5}=0;
%             else
%                 Model{i,5}=R;
%             end
%         end
%         
%         
%     end
 
end

