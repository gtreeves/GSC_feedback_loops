function signal_features = characterize_signal(X,t)

%here X can be a multi-channel
response = stepinfo(X,t);

for i=1:size(X,2)

    [pks,locs,widths,proms] = findpeaks(X(:,i),t);
    response(i).pks = pks;
    response(i).locs = locs;
    response(i).widths = widths;
    response(i).proms = proms;
    [a,b] = max(X(:,i));
    [c,d] = min(X(:,i));
    response(i).OS_calc = a - X(end,i);
    response(i).US_calc = c - X(end,i);
    response(i).OS_calc_time = t(b);
    response(i).US_calc_time = t(d);
end

signal_features = response;

end