function V = splitOnEpochStart(epoch)
    V = epoch.protocolSettings('epoch:startTime');
    V = V.getTime();
end