function V = splitOnBlockStart(epoch)
    V = epoch.protocolSettings('epochBlock:startTime');
    if ~isempty(V)
        V = V.getTime();
    end
end