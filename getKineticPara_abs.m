function [Kvec,realPro_]=getKineticPara_abs(Info_m,Asub,Aenz,Kvec_raw,Pro_,realMetab,realEconc)
%extracting the absolute K vector with input realEconc and realMetab

%-------------------------------
%Yikun Tan, UCLA
%Last revision Jan 14, 2010
%-------------------------------

nMetab_=size(Asub,2);
nPro_=size(Aenz,2);
totalRxns=size(Info_m,1);  %inclusing the inhibition reactions
nvuni_=Info_m(totalRxns,4);
nin_=length(find(Info_m(:,7)==2));
nout_=length(find(Info_m(:,7)==3));
nbm_=length(find(Info_m(:,7)==7));
nInterRxns=length(find(Info_m(:,7)==1|Info_m(:,7)==9));
nRevRxns=nInterRxns+nin_;
nRevEnz=max(Info_m(1:nRevRxns,2));
ntotal_Enz=size(unique(Info_m(:,2)),1);

for i=1:ntotal_Enz
    vrxns=find(Info_m(:,2)==i);
    for kk=1:length(vrxns)
        ind_tmp=[Info_m(vrxns(kk),3):1:Info_m(vrxns(kk),4)]; %row vector
        Kvec(ind_tmp,1)=Kvec_raw(ind_tmp)./repmat(realEconc(i),length(ind_tmp),1);
    end
end
subterm_real=exp([Asub Aenz]*log([realMetab;ones(nPro_,1)]));
Kvec=Kvec./subterm_real;


%calculate ref absolute enzyme distribution
realPro_=zeros(nPro_,1);
for i=1:nRevEnz
    vrxns=[];
    vrxns=find(Info_m(:,2)==i);
    for kk=1:length(vrxns)
        ind_tmp=[];
        ind_tmp=[Info_m(vrxns(kk),5):1:Info_m(vrxns(kk),6)]; %row vector
        realPro_(ind_tmp)=Pro_(ind_tmp).*repmat(realEconc(i),length(ind_tmp),1);
    end
    realPro_(i)=Pro_(i)*realEconc(i);
end
TransOut=find(Info_m(:,7)==3);
BMrxn=find(Info_m(:,7)==7);

if ~isempty(BMrxn)
 realPro_(Info_m(TransOut,5))=Pro_(Info_m(TransOut,5)).*realEconc(nRevEnz+1:1:end-1);
else
    %Change to realEconc(nRevEnz+2:1:end) for 2 enzymes catalyzing same reaction (14/29).
 realPro_(Info_m(TransOut,5))=Pro_(Info_m(TransOut,5)).*realEconc(nRevEnz+29:1:end);
end

