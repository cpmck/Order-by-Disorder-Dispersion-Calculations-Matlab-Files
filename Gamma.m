function val = Gamma( vi, vj, vs )
%Gamma() calculates Tr[M] 
%   calculate Gamma for transition between levels vi and vj of site s

%evaluate from base script
Jx = evalin('base','Jx');
Jy = evalin('base','Jy');
Jz = evalin('base','Jz');

%evaluate occupation probabilites from caller function
nA = evalin('caller','nA');

%choice of vs (1-4) returns systemA or systemB 
%1 and 3 for site A, 2 and 4 for site B
unitCell = zeros(length(Jx),length(Jx),4);

%evaluate energy eigenvectos from caller function
ssystemACalc = evalin('caller','ssystemACalc');
ssystemBCalc = evalin('caller','ssystemBCalc');

%contruct 9X9 matrices enumerated 1-4 for sites
unitCell(:,:,1) = ssystemACalc(:,2:end);
unitCell(:,:,2) = ssystemBCalc(:,2:end);
unitCell(:,:,3) = ssystemACalc(:,2:end);
unitCell(:,:,4) = ssystemBCalc(:,2:end);

%evaluate mean field moments from caller function
mm = zeros(1,3,4);
momCalc = evalin('caller','momCalc');

%construct 1X3 vectors enumerated by site vs
mm(1,:,1) = momCalc(1,:);
mm(1,:,2) = momCalc(2,:);
mm(1,:,3) = momCalc(1,:);
mm(1,:,4) = momCalc(2,:);

%select transition eignevector vi to vj on site vs
ei = unitCell(vi,:,vs)';
ef = unitCell(vj,:,vs)';

%select moment of site vs
mom = mm(1,:,vs);

%calculate gamma = tr[M] equation...
val1 = abs(ef'*(Jx- mom(1,1)*eye(length(Jx)))*ei).^2; 
val2 = abs(ef'*(Jy- mom(1,2)*eye(length(Jx)))*ei).^2;
val3 = abs(ef'*(Jz- mom(1,3)*eye(length(Jx)))*ei).^2;

%output Gamma with boltzmann weight
val = (val1+val2+val3)*(nA(vi,1)-nA(vj,1));

end

