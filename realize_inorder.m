function ew = realize_inorder(Db,Ds)
    ew = eig(Ds,Db); [~,ewidx] = sort(abs(ew),"descend"); ew = ew(ewidx);
end
