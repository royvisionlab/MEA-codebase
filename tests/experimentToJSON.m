function experimentToJSON(E, outdir, experimentDate)

% experimentDate = '20220420C';
% 
% E = experiment;

% Traverse the tree.
for pp = 1 : length(E.protocol)
    for gg = 1 : length(E.protocol(pp).group)
        if isfield(E.protocol(pp).group(gg),'block')
            for bb = 1 : length(E.protocol(pp).group(gg).block)
                for ee = 1 : length(E.protocol(pp).group(gg).block(bb).epoch)
                    k = keys(E.protocol(pp).group(gg).block(bb).epoch(ee).parameters);
                    for kk = 1 : length(k)
                        foo = E.protocol(pp).group(gg).block(bb).epoch(ee).parameters(k{kk});
%                         c = class(E.protocol(pp).group(gg).block(bb).epoch(ee).parameters(k{kk}));
                        if isa(foo, 'java.util.ArrayList')
                            tmp = convertArraylist(foo);
                            E.protocol(pp).group(gg).block(bb).epoch(ee).parameters(k{kk}) = tmp;
                        elseif isa(foo, 'java.util.Date')
                            E.protocol(pp).group(gg).block(bb).epoch(ee).parameters(k{kk}) = char(foo.toString);
                        end
                    end
                end
            end
        end
    end
end

txt = jsonencode(E);

% Save the text to a json file.
fid = fopen([outdir,experimentDate,'.json'], 'w');
fprintf(fid, '%s', txt);
fclose(fid);
end

function mArray = convertArraylist(foo)
    iterator = foo.iterator;
    mArray = zeros(1,foo.size);
    count=0;
    while (iterator.hasNext)
        count=count+1;
        mArray(count) = double(iterator.next);
    end
end