function indicator = checkstop(old,new,dt)
% indicate whether we should perform further iteraions or stop

    ind = find(abs(new)<=.5);
    M = length(ind);
    Q = sum(abs(new(ind)-old(ind)))./M;
    if Q<=dt*.18^2
        indicator = 1;
    else
        indicator = 0;
    end
end

