function merged=MergeColumns(a,b)
    if size(a,1)==size(b,1) && ...
            (size(a,2)==size(b,2) || size(a,2)==size(b,2)+1)
        merged=zeros(size(b,1),size(a,2)+size(b,2),class(a));
        merged(:,1:2:end)=a;
        merged(:,2:2:end)=b;
    end
end        
