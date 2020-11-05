function Rref=samplingR(Info_m,FElow,FEhigh,method_)

%Sample reversibilities under thermodynamic constraints

%-------------------------------
%Linh Tran, UCLA
%Matthew Rizk, UCLA
%Last revision March 27, 2009
%-------------------------------

totalRxns=size(Info_m,1);
nvuni_=Info_m(totalRxns,4);
nin_=length(find(Info_m(:,7)==2));
nout_=length(find(Info_m(:,7)==3));
nbm_=length(find(Info_m(:,7)==7));
nNetRxns=length(FElow);
nInterEnz=nNetRxns-nin_-nout_;
nRevRxns=nNetRxns-nout_;

tmpE=[FElow,FEhigh];
tmpE=abs(tmpE);

bound_(:,1)=min(tmpE,[],2);
bound_(:,2)=max(tmpE,[],2);
Rref=zeros(0.5*(nvuni_+nout_+nbm_),1);
rcount=1;
for i=1:nInterEnz+1:1:nRevRxns
    nsteps_=0.5*(Info_m(i,4)-Info_m(i,3)+1);
    d1=repmat(bound_(i,1),nsteps_,1);
    d2=repmat(bound_(i,2),nsteps_,1);
    y=1/nsteps_.*(-d1+(d1-d2).*(1-betarnd(4.5,1.5,nsteps_,1)));
    Rref(rcount:rcount+nsteps_-1,1)=exp(y);
    rcount=rcount+nsteps_;
end

rcount=rcount+nout_+nbm_;  %skip the outfluxes which have R=0;

if nNetRxns <= totalRxns
    Rref(rcount:end)=1;  %Reversibility of regulation reach equilibrium R=1;
end
nRegulationBinding=length(find(Info_m(:,7)==5|Info_m(:,7)==6|Info_m(:,7)==8));
nActivatedRxns=length(find(Info_m(:,7)==10));
for i=nNetRxns+nbm_+1:1:totalRxns
    if Info_m(i,7)==10
         nsteps_=0.5*(Info_m(i,4)-Info_m(i,3)+1);
         count_temp=0.5*(Info_m(i,3)+1+nout_+nbm_); %round(0.5*(Info_m(i,3)))+round(0.5*(nout_+nbm_));
         d1=repmat(bound_(Info_m(i,1),1),nsteps_,1);
         d2=repmat(bound_(Info_m(i,1),2),nsteps_,1);
         y=1/nsteps_.*(-d1+(d1-d2).*(1-betarnd(4.5,1.5,nsteps_,1)));
         Rref(count_temp:count_temp+nsteps_-1,1)=exp(y);
    end
end
        
    
