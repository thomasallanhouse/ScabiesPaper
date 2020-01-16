parfor i = 1:100
    [tim(:,i),y(:,:,i)] = SEIR_Gillespie(10,[48 0 2 0],0.4,0.1,0.2,0.5);
end