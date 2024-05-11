function archive = updateArchive(archive, pop, funvalue)


if archive.NP == 0, return; end

if size(pop, 1) ~= size(funvalue,1), error('check it'); end

popAll = [archive.pop;pop];

funvalues = [archive.funvalues; funvalue];
[dummy IX] = unique(popAll,'rows');
if length(IX) < size(popAll,1)
    popAll = popAll(IX,:);
    funvalues = funvalues(IX,:);
end

if size(popAll, 1) <= archive.NP   % add all new individuals
  archive.pop = popAll;
  archive.funvalues = funvalues;
else                % randomly remove some solutions
  rndpos = randperm(size(popAll, 1)); % equivelent to "randperm";
  rndpos = rndpos(1 : archive.NP);
  
  archive.pop = popAll  (rndpos, :);
  archive.funvalues = funvalues(rndpos, :);
end

