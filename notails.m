% This function removes the tails so that the plotting color makes more
% sense
function [largest,smallest,newdat] = notails(data,pctn)

    sorted = sort(data(:),'descend');
    largest = prctile(sorted,pctn);
    smallest = prctile(sorted,100-pctn);
    newdat = data(data<largest&data>smallest);
%     newdat = data;
end