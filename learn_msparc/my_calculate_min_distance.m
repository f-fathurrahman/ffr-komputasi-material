function nndis = my_calculate_min_distance(S)
if S.n_atm == 1
    pairs = [1, 1];
else 
    pairs = nchoosek(1:S.n_atm,2);
end
p1 = S.Atoms(pairs(:,1),:);
p2 = S.Atoms(pairs(:,2),:);

dd = calculateDistance(p1(:,1),p1(:,2),p1(:,3),p2(:,1),p2(:,2),p2(:,3),S);
nndis = min(dd);

end % function
