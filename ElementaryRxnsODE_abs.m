function [vuni,vnet,conc,rate_,t,initial_,LHS_]=ElementaryRxnsODE_abs(Kvec,Econc,Tconc,realEconc,Sp,Ae,Info_m,Vcof_index,refMetab,refPro,maxtime,step)

%ODE parameter
if isempty(maxtime)
    maxtime=1000;
    tmax = maxtime;
else
    tmax = maxtime;
end
if isempty(step)
    step_=maxtime/10;
else
    step_=step;
end

nin_=length(find(Info_m(:,7)==2));
nout_=length(find(Info_m(:,7)==3));
indInterEnzID=find(Info_m(:,7)==1 | Info_m(:,7)==9);
nInterEnz=length(unique(Info_m(indInterEnzID,2)));  %might be different from number of rxns
nInterRxn=length(find(Info_m(:,7)==1 | Info_m(:,7)==9));
nRevRxns=nInterRxn+nin_;
nNetRxns=nRevRxns+nout_;
totalRxns=size(Info_m,1);
if ~isempty(Vcof_index)
    temp_ind=Vcof_index(1);
else
    temp_ind=[];
end

%Find how namy complexes and free enzyme
if nNetRxns<totalRxns
    nEnzComp_=Info_m(totalRxns,6);
else
    nEnzComp_=Info_m(nRevRxns,6);
end
nRevEnz=max(Info_m(1:nRevRxns,2));  %no transporters for Vout
nMetab=size(Sp,1)-nEnzComp_;
nvuni_=size(Sp,2);

%%  Initialize the initial condition for ODE by setting
%%(1)all reversible enzymes and their intermediates: substrates
%%(2)transporters for outfluxes: included in k values
xo_=repmat(1,nin_,1);  %concentrations of the pathway substrates, which are constant and included into k values
EnzComp_=zeros(nEnzComp_,1); 
Ke=Kvec;

%initialize the metabolite concentrations
if isempty(refMetab)  
    refMetab=repmat(1,nMetab,1);
end
iniMetab=refMetab;

%initialize the enzyme concentrations
if isempty(refPro)
    flag_=1;    
else
    flag_=0;
    EnzComp_=refPro(1:nEnzComp_,1);
end

%Find the conservation pattern of cofactors
Nl=null(Sp','r');
nCofactors=length(Vcof_index);
allSel=[];
for jj=1:nCofactors
    seltmp=find(Nl(Vcof_index(jj),:)==-1);
    allSel=[allSel;seltmp'];
end
if length(allSel)>nCofactors
    allSel_new=unique(allSel);
    tmpCofac=Nl(Vcof_index,allSel_new);
    Xalpha=tmpCofac\(-eye(nCofactors));
    tmpCofac=[];
    tmpCofac=Nl(:,allSel_new)*Xalpha;
    Nl(:,allSel_new)=tmpCofac;
end
selCol=[];
for  jj=1:nCofactors
    selCol(jj)=find(Nl(Vcof_index(jj),:)==-1);
end
blockPatt=Nl(nMetab+1:end,selCol);
bookKeeping_=zeros(nRevEnz,nCofactors);


%Initialize internal enzymes
Econc=Econc.*realEconc(1:nInterEnz); %convert fold change to realEconc change
for k=1:nInterEnz
    vrxns=find(Info_m(:,2)==k); %find all including inhibition
    ind_tmp=k;
    for kk=1:length(vrxns)
        if Info_m(vrxns(kk),7)~=3  %not an Vout
            ind_tmp=[ind_tmp,Info_m(vrxns(kk),5):1:Info_m(vrxns(kk),6)]; %row vector
        end
    end
    
    if flag_
        tmp_nor=rand(length(ind_tmp),1);
        tmp_nor=tmp_nor/sum(tmp_nor);
        EnzComp_(ind_tmp,1)=tmp_nor;        
        tmp_E=1;
    else
        tmp_E=sum(EnzComp_(ind_tmp,1));
    end
    delE_=Econc(k,1)-tmp_E;  %free enzyme
    if delE_>0  %overexpress
        EnzComp_(k,1)=EnzComp_(k,1)+delE_;  %increase the free enzyme
    else   %underexpress
        EnzComp_(ind_tmp,1)=EnzComp_(ind_tmp,1)*Econc(k,1)/tmp_E;  %decrease each form by factor (total new/total old)
        bookKeeping_(k,:)=(1-tmp_E/Econc(k,1))*EnzComp_(ind_tmp,1)'*blockPatt(ind_tmp,:);
        for cc=1:nCofactors
            if bookKeeping_(k,cc)
                bookKeeping_(k,cc)=bookKeeping_(k,cc)-delE_;
            end
        end
    end
end

%Initialize transporters of influxes
Tconc=Tconc.*realEconc(nInterEnz+1:1:nInterEnz+nin_+nout_); %convert fold change to realEconc change
for k=nInterEnz+1:1:nInterEnz+nin_ 
    vrxns=find(Info_m(:,2)==k); %find all including inhibition
    rxnID=unique(Info_m(vrxns,1));
    if length(rxnID)>1
        error('\nPlease use separate transporter for different external metabolite');
    end
    ind_tmp=k;
    for kk=1:length(vrxns)
        if Info_m(vrxns(kk),7)~=3  %not an Vout
            ind_tmp=[ind_tmp,Info_m(vrxns(kk),5):1:Info_m(vrxns(kk),6)]; %row vector
        end
    end
    if flag_
        tmp_nor=rand(length(ind_tmp),1);
        tmp_nor=tmp_nor/sum(tmp_nor);
        EnzComp_(ind_tmp,1)=tmp_nor;        
        tmp_E=1;
    else
        tmp_E=sum(EnzComp_(ind_tmp,1,1));
    end
    delE_=Tconc(k-nInterEnz,1)-tmp_E;
    if delE_>0  %overexpress
        EnzComp_(k,1)=EnzComp_(k,1)+delE_;  %increase the free enzyme
    else   %underexpress
        EnzComp_(ind_tmp,1)=EnzComp_(ind_tmp,1)*Tconc(k-nInterEnz,1)/tmp_E;  %decrease each form by factor (total new/total old)
        bookKeeping_(k,:)=(1-tmp_E/Tconc(k-nInterEnz,1))*EnzComp_(ind_tmp,1)'*blockPatt(ind_tmp,:);
        for cc=1:nCofactors
            if bookKeeping_(k,cc)
                bookKeeping_(k,cc)=bookKeeping_(k,cc)-delE_;
            end
        end
    end
    Ke(Info_m(rxnID,3),1)=xo_(k-nInterEnz)*Ke(Info_m(rxnID,3),1);
end

%Check the input and initial outflux transporter
Vout_index=find(Info_m(:,7)==3);
if length(Vout_index)~=length(Tconc)-nin_
    error('\nPlease check the number of transporters');
end
Ke(Info_m(Vout_index,3),1)=Ke(Info_m(Vout_index,3),1).*Tconc(nin_+1:end);

if ~isempty(Vcof_index)
    iniMetab(Vcof_index)=refMetab(Vcof_index)+sum(bookKeeping_)';
end

%Initial concentrations
initial_=[iniMetab' EnzComp_']; 
[t,conc] = ode15s(@metabODE2,[0:step_:tmax],initial_,[],Ke,Ae,Sp,temp_ind);
rate_=rateCalc(conc,Ke,Ae);
vuni=rate_;
LHS_=norm(Sp*vuni(:,end)); %check the steady state
ind_for1=Info_m(1:nRevRxns,3);
ind_rev1=ind_for1+1;
for n=1:size(vuni,2)
    vnet(n,:)=[(vuni(ind_for1,n)-vuni(ind_rev1,n));vuni(Info_m(nRevRxns+1,3):Info_m(nRevRxns+nout_,3),n)]';
    %%%modify vnet to account for total flux through activated rxns
    %%%for Net15_3
    vnet(n,8)=[vuni(43,n)-vuni(44,n)+vuni(297,n)-vuni(298,n)];
    vnet(n,9)=[vuni(49,n)-vuni(50,n)+vuni(305,n)-vuni(306,n)];
    
    %%%for Net15_4
    %vnet(n,8)=[vuni(43,n)-vuni(44,n)+vuni(305,n)-vuni(306,n)];
    %vnet(n,9)=[vuni(49,n)-vuni(50,n)+vuni(313,n)-vuni(314,n)];
    
    %%%for Net15_6
    %vnet(n,40)=[vuni(235,n)-vuni(236,n)+vuni(305,n)-vuni(306,n)+vuni(313,n)-vuni(314,n)+vuni(321,n)-vuni(322,n)+vuni(329,n)-vuni(330,n)+vuni(337,n)-vuni(338,n)+vuni(345,n)-vuni(346,n)+vuni(353,n)-vuni(354,n)+vuni(361,n)-vuni(362,n)];
    
    %%%for Net15_7
    %vnet(n,40)=[vuni(235,n)-vuni(236,n)+vuni(305,n)-vuni(306,n)+vuni(313,n)-vuni(314,n)+vuni(321,n)-vuni(322,n)+vuni(329,n)-vuni(330,n)+vuni(337,n)-vuni(338,n)+vuni(345,n)-vuni(346,n)+vuni(353,n)-vuni(354,n)+vuni(361,n)-vuni(362,n)+vuni(369,n)-vuni(370,n)+vuni(377,n)-vuni(378,n)+vuni(385,n)-vuni(386,n)+vuni(393,n)-vuni(394,n)+vuni(401,n)-vuni(402,n)+vuni(409,n)-vuni(410,n)+vuni(417,n)-vuni(418,n)+vuni(425,n)-vuni(426,n)];
    
    %%%for Net15_8 to Net15_14
    %vnet(n,40)=[vuni(235,n)-vuni(236,n)+vuni(305,n)-vuni(306,n)+vuni(313,n)-vuni(314,n)+vuni(321,n)-vuni(322,n)+vuni(329,n)-vuni(330,n)+vuni(337,n)-vuni(338,n)+vuni(345,n)-vuni(346,n)+vuni(353,n)-vuni(354,n)+vuni(361,n)-vuni(362,n)];
    
    %%%for Net15_5
    %vnet(n,40)=[vuni(235,n)-vuni(236,n)+vuni(305,n)-vuni(306,n)+vuni(313,n)-vuni(314,n)+vuni(321,n)-vuni(322,n)+vuni(329,n)-vuni(330,n)];
    
    %%%for Net15_7_9
    %vnet(n,40)=[vuni(235,n)-vuni(236,n)+vuni(297,n)-vuni(298,n)+vuni(305,n)-vuni(306,n)+vuni(313,n)-vuni(314,n)+vuni(321,n)-vuni(322,n)+vuni(329,n)-vuni(330,n)+vuni(337,n)-vuni(338,n)+vuni(345,n)-vuni(346,n)+vuni(353,n)-vuni(354,n)+vuni(361,n)-vuni(362,n)+vuni(369,n)-vuni(370,n)+vuni(377,n)-vuni(378,n)+vuni(385,n)-vuni(386,n)+vuni(393,n)-vuni(394,n)+vuni(401,n)-vuni(402,n)+vuni(409,n)-vuni(410,n)+vuni(417,n)-vuni(418,n)];
    
    %%%for Net15_7_10
    %vnet(n,8)=[vuni(43,n)-vuni(44,n)+vuni(297,n)-vuni(298,n)];
    %vnet(n,9)=[vuni(49,n)-vuni(50,n)+vuni(305,n)-vuni(306,n)];
    %vnet(n,40)=[vuni(235,n)-vuni(236,n)+vuni(313,n)-vuni(314,n)+vuni(321,n)-vuni(322,n)+vuni(329,n)-vuni(330,n)+vuni(337,n)-vuni(338,n)+vuni(345,n)-vuni(346,n)+vuni(353,n)-vuni(354,n)+vuni(361,n)-vuni(362,n)+vuni(369,n)-vuni(370,n)+vuni(377,n)-vuni(378,n)+vuni(385,n)-vuni(386,n)+vuni(393,n)-vuni(394,n)+vuni(401,n)-vuni(402,n)+vuni(409,n)-vuni(410,n)+vuni(417,n)-vuni(418,n)+vuni(425,n)-vuni(426,n)+vuni(433,n)-vuni(434,n)];
    
    %%%for Net15_14_11
    %vnet(n,6)=[vuni(31,n)-vuni(32,n)+vuni(297,n)-vuni(298,n)];
    %vnet(n,7)=[vuni(37,n)-vuni(38,n)+vuni(305,n)-vuni(306,n)];
    
    %%%for Net16_1
    %vnet(n,8)=[vuni(43,n)-vuni(44,n)+vuni(361,n)-vuni(362,n)];
    %vnet(n,9)=[vuni(49,n)-vuni(50,n)+vuni(369,n)-vuni(370,n)];
    
    %%%for Net16_2
    %vnet(n,8)=[vuni(43,n)-vuni(44,n)+vuni(361,n)-vuni(362,n)];
    %vnet(n,9)=[vuni(49,n)-vuni(50,n)+vuni(369,n)-vuni(370,n)];
    
    %%%for Net16_3 & Net16_4
    %vnet(n,8)=[vuni(43,n)-vuni(44,n)+vuni(297,n)-vuni(298,n)];
    %vnet(n,9)=[vuni(49,n)-vuni(50,n)+vuni(305,n)-vuni(306,n)];
    
    %%%for Net16_5
    %vnet(n,8)=[vuni(43,n)-vuni(44,n)+vuni(329,n)-vuni(330,n)];
    %vnet(n,9)=[vuni(49,n)-vuni(50,n)+vuni(337,n)-vuni(338,n)];
end
%-----------------------------------------------------------

function dx = metabODE2(t,x,Ke,Ae,Sp,temp_ind)

[nrxns,nsubs]=size(Ae);
dx = zeros(nsubs,1);
V=zeros(nrxns,1);
log_x=log(x);
tmp=sign(Ae)*log_x;
V=Ke.*exp(tmp);
dx=Sp*V;
% d_atp=-4163/1000/(7/20+5/4*t+1/20*t^2)+4163/1000*t/(7/20+5/4*t+1/20*t^2)^2*(5/4+1/10*t);
% d_adp=-519/2000*(2731/1000)^(-3/20*t)*log(2731/1000)*(3/25*t+986900807943461/4611686018427387904*t^3)+173/100*(2731/1000)^(-3/20*t)*(3/25+2960702423830383/4611686018427387904*t^2);
% d_nadh=-1353/10000000*(2371/1000)^(-123/1000*t)*log(2371/1000)*(211/250*t+13/125*t^3)+11/10000*(2371/1000)^(-123/1000*t)*(211/250+39/125*t^2);
% d_nad=-57159/1000000*(273/100)^(-87/2000*t-171/500)*log(273/100)-(273/100)^(-109/5000*t-171/1000)/(8481/1000+t)+109/5000*(t+7871/1000)*(273/100)^(-109/5000*t-171/1000)*log(273/100)/(8481/1000+t)+(t+7871/1000)*(273/100)^(-109/5000*t-171/1000)/(8481/1000+t)^2;
% d_nadph=-2407/15625*(1359/500)^(-58/125*t)*log(1359/500)*(83/5000*t^(79/50)+6124319032471571/36893488147419103232*t^(473/100)+4371475763726499/38685626227668133590597632*t^(789/100)+673439381371247/4951760157141521099596496896*t^11+4989472762498311/40564819207303340847894502572032*t^(71/5))+83/250*(1359/500)^(-58/125*t)*(6557/250000*t^(29/50)+2896802902359053083/3689348814741910323200*t^(373/100)+3449094377580207711/3868562622766813359059763200*t^(689/100)+7407833195083717/4951760157141521099596496896*t^10+354252566137380081/202824096036516704239472512860160*t^(66/5));
% d_nadp=14167099448608935/2361183241434822606848*t^2+1/2500*t-1/100;
% dx(temp_ind)=d_atp;
% dx(temp_ind+1)=d_adp;
% dx(temp_ind+2)=d_nadh;
% dx(temp_ind+3)=d_nad;
% dx(temp_ind+4)=d_nadph;
% dx(temp_ind+5)=d_nadp;
%----------------------------------------------
function rate_=rateCalc(xconc,Ke,Ae)

ndata=size(xconc,1);
[nrxns,nsubs]=size(Ae);
rate_=zeros(nrxns,ndata);
log_x=log(xconc)';
tmp=sign(Ae)*log_x;
for k=1:ndata
    rate_(:,k)=Ke.*exp(tmp(:,k));
end


