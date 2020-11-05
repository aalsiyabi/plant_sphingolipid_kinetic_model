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