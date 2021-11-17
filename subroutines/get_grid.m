function [grid] = get_grid(gridvals,no_traits)
if isrow(gridvals)
    gridvals=gridvals';
end
gridvals = unique(abs(gridvals(gridvals~=0)));
gridvals=[-gridvals; gridvals; 0];
gridsize=length(gridvals);

no_combos=gridsize^(no_traits);
grid=zeros(no_traits,no_combos);
for kk=1:no_traits
    row=repmat(gridvals,gridsize^(no_traits-kk),gridsize^(kk-1))';
    grid(kk,:)=row(:);
end

grid=grid(:,grid(1,:)>=0);
grid=grid(:,~all(grid==0,1));
end