function [Spsub,Spenz]=getExpandedS(S,Sreg,Sbm,Vin_index,Vout_index,Info_m)

%Converting the stoichiometric coefficient matrix of overall reactions (S)
%and the stoichiometric biomass matrix (Sbm) into the corresponding
%stoichiometric coefficient matrix of elementary reactions.  
%The transporters corresponding to the input and output must be last columns
%in S. The influx is modeled as reversible, while the outflux(es) is(are) 
%irreversible.
%
%-------------------------------
%Linh Tran, UCLA
%Matthew Rizk, UCLA
%Yikun Tan, UCLA
%Last revision Jan 14, 2010
%-------------------------------

[nMetab,nNetRxns]=size(S);


nin_=length(Vin_index);
nout_=length(Vout_index);
nInterRxns=nNetRxns-nin_-nout_;
nRevRxns=nInterRxns+nin_;
totalRevEnz=length(unique(Info_m(1:nRevRxns,2)));

vtmp=find(S>0);
Spro=zeros(nMetab,nNetRxns);
Spro(vtmp)=S(vtmp);

vtmp=find(S<0);
Ssub=zeros(nMetab,nNetRxns);
Ssub(vtmp)=S(vtmp);

vnsub=sum(Ssub);
vnpro=sum(Spro);
rcount=1+totalRevEnz;
ccount=1;
Spenz=[];
Spsub=[];

vtmp=find(Sbm>0);
Sbmpro=zeros(size(Sbm));
Sbmpro(vtmp)=Sbm(vtmp);

vtmp=find(Sbm<0);
Sbmsub=zeros(size(Sbm));
Sbmsub(vtmp)=Sbm(vtmp);

vnbmsub=sum(Sbmsub);
vnbmpro=sum(Sbmpro);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INTERNAL REACTIONS

for k=1:nInterRxns
    if vnsub(k)==-2 && vnpro(k)==2
    %two substrates and two products: 4 steps, 4 enzyme complexes, and 8 Vuni 
        blockSsub=zeros(nMetab,8);
        subID=getExpandedIndices(Ssub(:,k));
        blockSsub(subID(1),1)=-1;
        blockSsub(subID(2),3)=-1;
        blockSsub(subID(1),2)=1;
        blockSsub(subID(2),4)=1;            
        proID=getExpandedIndices(Spro(:,k)); 
        blockSsub(proID(1),6)=-1;
        blockSsub(proID(2),8)=-1;
        blockSsub(proID(1),5)=1;
        blockSsub(proID(2),7)=1;               
        blockSenz=zeros(4,8);
        blockSenz(1,[1 8])=-1;
        blockSenz(2,2:3)=-1; 
        blockSenz(3,4:5)=-1;
        blockSenz(4,6:7)=-1;
        blockSenz(1,[2 7])=1;
        blockSenz(2,[1 4])=1; 
        blockSenz(3,[3 6])=1;
        blockSenz(4,[5 8])=1;
    else
        blockSenz=zeros(3,6);
        blockSenz(1,[1 6])=-1;
        blockSenz(2,2:3)=-1; 
        blockSenz(3,4:5)=-1;
        blockSenz(1,[2 5])=1;
        blockSenz(2,[1 4])=1; 
        blockSenz(3,[3 6])=1;
        blockSsub=zeros(nMetab,6);
        subID=getExpandedIndices(Ssub(:,k));
        if vnsub(k)==-2 & vnpro(k)==1
        %2 substrates and 1 products: 3 steps,3 enzyme complexes, and 6 Vuni    
            blockSsub(subID(1),1)=-1;
            blockSsub(subID(2),3)=-1;
            blockSsub(subID(1),2)=1;
            blockSsub(subID(2),4)=1;
            proID=getExpandedIndices(Spro(:,k));
            blockSsub(proID,6)=-1;
            blockSsub(proID,5)=1;
        else  %1 substrates        
            blockSsub(subID,1)=-1;
            blockSsub(subID,2)=1;
            proID=getExpandedIndices(Spro(:,k));
            if vnpro(k)==1  %1 product
                blockSsub(proID,6)=-1;
                blockSsub(proID,5)=1;
            else             %2 products
                blockSsub(proID(1),4)=-1;
                blockSsub(proID(2),6)=-1;
                blockSsub(proID(1),3)=1;
                blockSsub(proID(2),5)=1;                
            end
        end
    end
    if vnsub(k)==-2 & vnpro(k)==3
        %two substrates and three products: 5 steps,5 enzyme complexes, and 10 Vuni
        blockSsub=zeros(nMetab,10);
        subID=getExpandedIndices(Ssub(:,k));
        blockSsub(subID(1),1)=-1;
        blockSsub(subID(2),3)=-1;
        blockSsub(subID(1),2)=1;
        blockSsub(subID(2),4)=1;
        proID=getExpandedIndices(Spro(:,k));
        blockSsub(proID(1),6)=-1;
        blockSsub(proID(2),8)=-1;
        blockSsub(proID(3),10)=-1;
        blockSsub(proID(1),5)=1;
        blockSsub(proID(2),7)=1;
        blockSsub(proID(3),9)=1;
        blockSenz=zeros(5,10);
        blockSenz(1,[1 10])=-1;
        blockSenz(2,2:3)=-1;
        blockSenz(3,4:5)=-1;
        blockSenz(4,6:7)=-1;
        blockSenz(5,8:9)=-1;
        blockSenz(1,[2 9])=1;
        blockSenz(2,[1 4])=1;
        blockSenz(3,[3 6])=1;
        blockSenz(4,[5 8])=1;
        blockSenz(5,[7 10])=1;
    end
    if vnsub(k)<-2 || vnpro(k)>3
        error('\nThe program limits to maximum of 2 reactants and 3 products per reaction\n');
    end
    Spsub=[Spsub,blockSsub];
    [nr,ne]=size(blockSenz);
    Spenz(Info_m(k,2),ccount:ccount+ne-1)=blockSenz(1,:);
    Spenz(rcount:rcount+nr-2,ccount:ccount+ne-1)=blockSenz(2:nr,:);
    rcount=rcount+nr-1;
    ccount=ccount+ne;
end             

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%UPTAKE REACTIONS

for k=1:nin_ 
    blockSsub=zeros(nMetab,6);
    subID=getExpandedIndices(Ssub(:,nInterRxns+k));
    blockSsub(subID,1)=-1;
    blockSsub(subID,2)=1;
    proID=getExpandedIndices(Spro(:,nInterRxns+k));
    if vnpro(nInterRxns+k)==1  %1 product
        blockSsub(proID,6)=-1;
        blockSsub(proID,5)=1;
    else             %2 products
        blockSsub(proID(1),4)=-1;
        blockSsub(proID(2),6)=-1;
        blockSsub(proID(1),3)=1;
        blockSsub(proID(2),5)=1;                
    end
    blockSenz=zeros(3,6);
    blockSenz(1,[1 6])=-1;
    blockSenz(2,2:3)=-1; 
    blockSenz(3,4:5)=-1;
    blockSenz(1,[2 5])=1;
    blockSenz(2,[1 4])=1; 
    blockSenz(3,[3 6])=1;
    Spsub=[Spsub,blockSsub];
    [nr,ne]=size(blockSenz);
    Spenz(Info_m(nInterRxns+k,2),ccount:ccount+ne-1)=blockSenz(1,:);
    Spenz(rcount:rcount+nr-2,ccount:ccount+ne-1)=blockSenz(2:nr,:); 
    rcount=rcount+nr-1;
    ccount=ccount+ne;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT REACTIONS

nTinEnz_=size(Spenz,1);
for k=1:nout_
    subID=find(Ssub(:,Vout_index(k)));
    blockSsub=zeros(nMetab,1); 
    blockSsub(subID,1)=-1;   
    blockSenz=zeros(nTinEnz_,1);
    Spsub=[Spsub,blockSsub];
    Spenz=[Spenz,blockSenz];  %add column corresponding to V, but not row for enz
    ccount=ccount+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BIOMASS REACTION
if ~isempty(Sbm)
    nTinEnz_=size(Spenz,1);
    subID=find(Sbmsub(:,1));
    proID=find(Sbmpro(:,1));
    blockSsub=zeros(nMetab,1);
    blockSsub(subID,1)=Sbm(subID,1);
    blockSsub(proID,1)=Sbm(proID,1);
    blockSenz=zeros(nTinEnz_,1);
    Spsub=[Spsub,blockSsub];
    Spenz=[Spenz,blockSenz];  %add column corresponding to V, but not row for enz
    ccount=ccount+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%REGULATION REACTIONS
if ~isempty(Sreg)
    [xreg,yreg]=find(Sreg<0);
    vtmp=find(Sreg<0);
    for j=1:length(vtmp)
        located_=[];
        nNormalEnz=rcount-1;
        if Sreg(vtmp(j))==-1
            %competitive inhibition:  1 step, 1 enzyme complex, and 2 Vuni
            blockSregsub=zeros(nMetab,2);
            blockSregsub(xreg(j),1)=-1;
            blockSregsub(xreg(j),2)=1;
            blockSregenz1=zeros(nNormalEnz,2);
            located_=Info_m(yreg(j),2);
            blockSregenz1(located_,:)=[-1 1];
            blockSregenz2=[1 -1];
        elseif Sreg(vtmp(j))==-2
            %uncompetitive inhibition:  1 step, 1 enzyme complex, and 2 Vuni
            blockSregsub=zeros(nMetab,2);
            blockSregsub(xreg(j),1)=-1;
            blockSregsub(xreg(j),2)=1;
            blockSregenz1=zeros(nNormalEnz,2);
            located_=Info_m(yreg(j),5);
            blockSregenz1(located_,:)=[-1 1];
            blockSregenz2=[1 -1];
        elseif Sreg(vtmp(j))==-3
            %mixed inhibition:  2 steps, 2 enzyme complexes, and 4 Vuni
            blockSregsub=zeros(nMetab,4);
            blockSregsub(xreg(j),[1 3])=-1;
            blockSregsub(xreg(j),[2 4])=1;
            blockSregenz1=zeros(nNormalEnz,4);
            located_(1)=Info_m(yreg(j),2);
            located_(2)=Info_m(yreg(j),5);
            blockSregenz1(located_(1),1:2)=[-1 1];
            blockSregenz1(located_(2),3:4)=[-1 1];
            blockSregenz2=zeros(2,4);
            blockSregenz2(1,1:2)=[1 -1];
            blockSregenz2(2,3:4)=[1 -1];
        elseif Sreg(vtmp(j))==-4               %Activated rxn and Activation rxn
            k=yreg(j);
            if vnsub(k)==-2 & vnpro(k)==2
                %two substrates and two products: 4 steps, 4 enzyme complexes, and 8 Vuni
                blockSregsub=zeros(nMetab,10);
                blockSregsub(xreg(j),1)=-1;
                blockSregsub(xreg(j),2)=1;
                subID=getExpandedIndices(Ssub(:,k));
                blockSregsub(subID(1),3)=-1;
                blockSregsub(subID(2),5)=-1;
                blockSregsub(subID(1),4)=1;
                blockSregsub(subID(2),6)=1;
                proID=getExpandedIndices(Spro(:,k));
                blockSregsub(proID(1),8)=-1;
                blockSregsub(proID(2),10)=-1;
                blockSregsub(proID(1),7)=1;
                blockSregsub(proID(2),9)=1;
                blockSregenz1=zeros(nNormalEnz,10);
                blockSregenz1(Info_m(yreg(j),2),1)=-1;
                blockSregenz1(Info_m(yreg(j),2),2)=1;
                blockSregenz2=zeros(4,10);
                blockSregenz2(1,[2 3 10])=-1;
                blockSregenz2(2,4:5)=-1;
                blockSregenz2(3,6:7)=-1;
                blockSregenz2(4,8:9)=-1;
                blockSregenz2(1,[1 4 9])=1;
                blockSregenz2(2,[3 6])=1;
                blockSregenz2(3,[5 8])=1;
                blockSregenz2(4,[7 10])=1;
            else
                blockSregenz1=zeros(nNormalEnz,8);
                blockSregenz1(Info_m(yreg(j),2),1)=-1;
                blockSregenz1(Info_m(yreg(j),2),2)=1;
                blockSregenz2=zeros(3,8);
                blockSregenz2(1,2:3)=-1;
                blockSregenz2(2,4:5)=-1;
                blockSregenz2(3,6:7)=-1;
                blockSregenz2(1,[1 4])=1;
                blockSregenz2(2,[3 6])=1;
                blockSregenz2(3,[5 8])=1;
                blockSregsub=zeros(nMetab,8);
                subID=getExpandedIndices(Ssub(:,k));
                if vnsub(k)==-2 & vnpro(k)==1
                    %2 substrates and 1 products: 3 steps,3 enzyme complexes, and 6 Vuni
                    blockSregsub(xreg(j),1)=-1;
                    blockSregsub(xreg(j),2)=1;
                    blockSregsub(subID(1),3)=-1;
                    blockSregsub(subID(2),5)=-1;
                    blockSregsub(subID(1),4)=1;
                    blockSregsub(subID(2),6)=1;
                    proID=getExpandedIndices(Spro(:,k));
                    blockSregsub(proID,8)=-1;
                    blockSregsub(proID,7)=1;
                else  %1 substrates
                    blockSregsub(xreg(j),1)=-1;
                    blockSregsub(xreg(j),2)=1;
                    blockSregsub(subID,3)=-1;
                    blockSregsub(subID,4)=1;
                    proID=getExpandedIndices(Spro(:,k));
                    if vnpro(k)==1  %1 product
                        blockSregsub(proID,8)=-1;
                        blockSregsub(proID,7)=1;
                    else             %2 products
                        blockSregsub(proID(1),6)=-1;
                        blockSregsub(proID(2),8)=-1;
                        blockSregsub(proID(1),5)=1;
                        blockSregsub(proID(2),7)=1;
                    end
                end
            end
            if vnsub(k)==-2 & vnpro(k)==3
                %two substrates and three products: 5 steps,5 enzyme complexes, and 10 Vuni
                blockSregsub=zeros(nMetab,12);
                subID=getExpandedIndices(Ssub(:,k));
                blockSregsub(xreg(j),1)=-1;
                blockSregsub(xreg(j),2)=1;
                blockSregsub(subID(1),3)=-1;
                blockSregsub(subID(2),5)=-1;
                blockSregsub(subID(1),4)=1;
                blockSregsub(subID(2),6)=1;
                proID=getExpandedIndices(Spro(:,k));
                blockSregsub(proID(1),8)=-1;
                blockSregsub(proID(2),10)=-1;
                blockSregsub(proID(3),12)=-1;
                blockSregsub(proID(1),7)=1;
                blockSregsub(proID(2),9)=1;
                blockSregsub(proID(3),11)=1;
                blockSregenz1=zeros(nNormalEnz,12);
                blockSregenz1(Info_m(yreg(j),2),1)=-1;
                blockSregenz1(Info_m(yreg(j),2),2)=1;
                blockSregenz2=zeros(5,12);
                blockSenz(1,[2 3 12])=-1;
                blockSenz(2,4:5)=-1;
                blockSenz(3,6:7)=-1;
                blockSenz(4,8:9)=-1;
                blockSenz(5,10:11)=-1;
                blockSenz(1,[1 4 11])=1;
                blockSenz(2,[3 6])=1;
                blockSenz(3,[5 8])=1;
                blockSenz(4,[7 10])=1;
                blockSenz(5,[9 12])=1;
            end
            
        else
            error('\nThe program only accepts inhibition labeled with -1 or -2 or -3\n');
        end
        Spsub=[Spsub,blockSregsub];
        [nr1,ne1]=size(blockSregenz1);
        [nr2,ne2]=size(blockSregenz2);
        Spenz=[Spenz,blockSregenz1];
        Spenz(rcount:rcount+nr2-1,ccount:ccount+ne2-1)=blockSregenz2;
        rcount=rcount+nr2;
        ccount=ccount+ne2;
    end
end

%------------------------------------------------------
function v_ind=getExpandedIndices(v_input)

n=length(v_input);
v_input=abs(v_input);
v_ind=[];
for k=1:n
    if v_input(k)==1
        v_ind=[v_ind;k];
    elseif v_input(k)>1
        v_ind=[v_ind;repmat(k,v_input(k),1)];
    end
end