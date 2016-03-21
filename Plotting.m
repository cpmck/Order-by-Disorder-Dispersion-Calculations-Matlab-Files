dispersion = dispersion1;

xpoints = [-1:0.01:1]';

dispe = zeros(length(h),1,length(dispersion));
dispE = zeros(length(h),1);
for k = 1:length(dispersion)
    for i = 1:length(h)
   
        dispe(i,1,k) = dispersion(k,1,i);
    end
    
    dispE = horzcat(dispE,dispe(:,1,k));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dispi = zeros(length(h),1,length(dispersion));
dispI = zeros(length(h),1);
for k = 1:length(dispersion)
    for i = 1:length(h)
   
        dispi(i,1,k) = dispersion(k,2,i);
    end
    dispI = horzcat(dispI,dispi(:,1,k));
end


pinten = @(j,x) interp1(h, dispI(:,j),x,'spline','extrap');
plin = @(j,x) interp1(h, dispE(:,j),x,'spline','extrap');



sim = @(j,x,y) pinten(j,x).*exp(-(y-plin(j,x)).^2 * (1/0.2));


simulation = @(x,y) sim(1,x,y);

for j = 2:length(dispersion)
   
    simulation = @(x,y) simulation(x,y) + sim(j,x,y);
    
end


[X,Y] = meshgrid(-1:.01:1,9:.01:20);
Z = simulation(X,Y);

figure
mesh(X,Y,Z)

figure
contourf(X,Y,Z,8);




