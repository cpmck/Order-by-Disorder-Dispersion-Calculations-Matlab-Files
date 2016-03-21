function Jreturn = JOp(x)
%Return Jx, Jy or Jz = J(1),J(2) and J(3)

%row vector of m values
setVals = [-4:1:4];

%calculate values of raising and lowering operators
Jp = zeros(9,9);
Jm = zeros(9,9);
for i = 1:9 

    for j = 1:9
        if j-i == 1
            Jp(i,j) = sqrt(20 - (i-5)*(i-4));
        elseif i-j == 1
            Jm(i,j) = sqrt(20 - (i-5)*(i-6));
        end
    end
end

setVals = [-4:1:4];

if x == 3
Jz = diag(setVals);
Jreturn = Jz;
elseif x == 1
Jx = 1/2*(Jp+Jm);
Jreturn = Jx;
elseif x == 2
Jy = 1/(2*1i) * (Jp-Jm);
Jreturn = Jy;
end

end

