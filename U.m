function Ureturn = U( vi, vj, vs )
%Ureturn returns column vector U for 2 level system vi and vj of site s
%   This function calculates U Eqn....

%evalueate angular momentum operators from base script
Jx = evalin('base','Jx');
Jy = evalin('base','Jy');
Jz = evalin('base','Jz');

%choice of vs (1-4) returns systemA or systemB
%1 and 3 for A site, 2 and 4 for B site
unitCell = zeros(length(Jx),length(Jx),4);

%evaluate from caller function
ssystemACalc = evalin('caller','ssystemACalc');
ssystemBCalc = evalin('caller','ssystemBCalc');

%input eigenvectors into 9X9 matrices enumerated 1-4 for the unit cell site
unitCell(:,:,1) = ssystemACalc(:,2:end);
unitCell(:,:,2) = ssystemBCalc(:,2:end);
unitCell(:,:,3) = ssystemACalc(:,2:end);
unitCell(:,:,4) = ssystemBCalc(:,2:end);

%import mean field moments and hold in mm enumerated 1-4 for unit cell site
mm = zeros(1,3,4);
momCalc = evalin('caller','momCalc');
mm(1,:,1) = momCalc(1,:);
mm(1,:,2) = momCalc(2,:);
mm(1,:,3) = momCalc(1,:);
mm(1,:,4) = momCalc(2,:);

%select initial and final state eignvectors 
ei = transpose(unitCell(vi,:,vs));
ef = transpose(unitCell(vj,:,vs));

%select moment 
mom = mm(1,:,vs);

%calculate norm squared values
val1 = abs(ef'*(Jx- mom(1,1)*eye(length(Jx)))*ei).^2; 
val2 = abs(ef'*(Jy- mom(1,2)*eye(length(Jx)))*ei).^2;
val3 = abs(ef'*(Jz- mom(1,3)*eye(length(Jx)))*ei).^2;

%construct denominator
denominator = sqrt(val1+val2+val3);

%calculate numerators
num1 = ei'*(Jx- mom(1,1)*eye(length(Jx)))*ef;
num2 = ei'*(Jy- mom(1,2)*eye(length(Jx)))*ef;
num3 = ei'*(Jz- mom(1,3)*eye(length(Jx)))*ef;

%return the column vector
Ureturn = [num1;num2;num3]/denominator;
end

