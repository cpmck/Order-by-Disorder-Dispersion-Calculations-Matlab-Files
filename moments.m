function returnj = moments( H,beta,phi,IAA,IAB )
%MOMENTS Calculate mean field moments
%   This function impliments a mean field self consistency calculation.
%   moments() uses the function newj() to calculate new moments at each
%   iteration. 


%set converegence limit
error = 10e-15;
%make an initial estimate of moments
j0 = [0,0;0.1,0.1;0,0];
j = j0;

%k is used to start the self consistency algorithm
k = [9,9;9,9;9,9];

%define newj() as a function of the moments only
fun = @(j) real(newj(j(:,1),j(:,2),H,beta,phi,IAA,IAB));

%Iterate until convergence is reached
while max( abs(j(:)-k(:)) ) > error    
    j = k;
    k = fun(j);  
end

%return moments
returnj = transpose(j);
end
