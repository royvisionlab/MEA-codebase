
wPCA  = results.ops.wPCA;
wTEMP  = results.ops.wTEMP;

tF = results.cProjPC;
Nchan = results.ops.Nchan;

ss = double(results.st3(:,1)) / results.ops.fs;

ktid = int32(results.st3(:,2));

[iC, mask, C2C] = getClosestChannels(results, results.ops.sigmaMask, min(results.ops.Nchan, 32));

[~,iW] = max(abs(results.dWU(results.ops.nt0min,:,:)),[],2);
iW = int32(squeeze(iW));

iC = gather(iC(:,iW));

uweigh = abs(results.U(:,:,1));
uweigh = uweigh ./ sum(uweigh,1);
ycup = sum(uweigh .* results.yc, 1);
xcup = sum(uweigh .* results.xc, 1);


dmin = median(diff(unique(results.yc)));

yunq = unique(results.yc);
mxc = zeros(numel(yunq), 1);
for j = 1:numel(yunq)
    xc = results.xc(results.yc==yunq(j));
    if numel(xc)>1
       mxc(j) = median(diff(sort(xc))); 
    end
end
dminx = max(5, median(mxc));

ycenter = (min(results.yc) + dmin-1):(2*dmin):(max(results.yc)+dmin+1);
xcenter = (min(results.xc) + dminx-1):(2*dminx):(max(results.xc)+dminx+1);
[xcenter, ycenter] = meshgrid(xcenter, ycenter);
xcenter = xcenter(:);
ycenter = ycenter(:);

Wpca = zeros(6, Nchan, 1000, 'single');
nst = numel(ktid);
hid = zeros(nst,1 , 'int32');

xy = zeros(nst, 2);

tic
for j = 1:numel(ycenter)
    y0 = ycenter(j);
    x0 = xcenter(j);    
    xchan = (abs(ycup - y0) < dmin) & (abs(xcup - x0) < dminx);
    
    itemp = find(xchan);
        
    tin = ismember(ktid, itemp);
    
    if sum(tin)<1
        continue;
    end
    
    pid = ktid(tin);
    data = tF(tin, :, :);
    
    ich = unique(iC(:, itemp));
%     ch_min = ich(1)-1;
%     ch_max = ich(end);
    
    nsp = size(data,1);
    dd = zeros(nsp, 3,  numel(ich),  'single');
%     dd = zeros(nsp, 6,  numel(ich),  'single');
    for k = 1:length(itemp)
        ix = pid==itemp(k);
        [~,ia,ib] = intersect(iC(:,itemp(k)), ich);
        dd(ix, :, ib) = data(ix,:,ia);
    end
    xy(tin, :) =  spike_position(dd, wPCA, wTEMP, results.xc(ich), results.yc(ich));
end

function xy = spike_position(dd, wPCA, wTEMP, yc, xc)
dT = gpuArray(dd);
[nspikes, nPC, nchan] = size(dT);
[~, imax] = max(max(dT.^2, [], 3), [], 2);
dBest = gpuArray.zeros(nspikes, nchan, 'single');
for j = 1:nPC
    iX = imax==j;
    dBest(iX, :) = dT(iX, j, :);
end

dBest = max(0, dBest);
dBest = dBest ./ sum(dBest,2);

ysp = dBest * yc;
xsp = dBest * xc;

xy = [xsp, ysp];

xy = gather(xy);
end