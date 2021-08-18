function wH = radialweightedHist(rinitial,rfinal,edges)

%The expontential can be added to change the profile

% wH = sum(exp(-1.5*(rinitial/max(rinitial)).^2).*(rinitial./rfinal).^1.*(rfinal>=edges(1:end-1)')...
%     .*(rfinal<edges(2:end)'),2)';

wH = sum((rinitial./(rfinal+1e-9)).^1.*(rfinal>=edges(1:end-1)')...
    .*(rfinal<edges(2:end)'),2)';

end