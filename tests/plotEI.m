
electrode = 348;



xmap = round(foo.electrode_map(:,1));
ymap = round(foo.electrode_map(:,2));

img = meshgrid(xmap,ymap);

tmp = foo.E.(['e',num2str(electrode)]);

for jj = 1 : size(tmp,2)
    im = img;
    im
end


