function gmd = gmatch(rew,cew)
    new = zeros(size(cew));
    % greedily match as many ew/ev as we can to reference
    for i=1:min(length(rew),length(cew))
        mdist = Inf; midx = -1;
        for j=1:length(cew)
            if norm(rew(i)-cew(j)) < mdist
                mdist = norm(rew(i)-cew(j)); midx = j;
            end
        end
        new(i) = cew(midx); cew(midx) = [];
    end
    gmd = norm(rew-new);
end