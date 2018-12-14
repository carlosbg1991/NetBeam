%% Apply Kriging and get prediction and confidence intervals for sampled data

condPoints = [1,2,3,4,5,6,7]';  % Sample Points
condVals = [5,6,7,6,4,5,4]';  % Sample values

lx = 10;  % Max value to predict
nx = 1e3;  % Maximum number of points to predict
dx = lx/nx;  % Distance between prediction points
uncondPoints = [dx/2:dx:lx]';  % Interpolated range

corFun = 'sexp';  % Kernel type: 'exp','sexp','poly','tri'
lowerTheta = 0;  % Lower correlation bound (range for samples to be corr)
upperTheta = 100;  % Upper correlation bound (range for samples to be corr)

% Maximum Likelihood estimation
[theta,mu,sigma,lval] = maxLfun(condVals,condPoints,corFun,lowerTheta,upperTheta);

% Interpolate over the unknown values using Kriging. Returns predictions
% and confidence bounds
[krige,CIupper,CIlower] = krigeIt(condPoints,condVals,uncondPoints,corFun,mu,sigma,theta);

%% Plotting section
figure(1)
scatter(condPoints,condVals)
hold on
plot(uncondPoints,krige)


fig2 = figure(2);
hold on
X = [uncondPoints; flipud(uncondPoints)];
Y = [CIlower; flipud(CIupper)];
h = fill(X',Y',[0.8,0.8,0.8]);
set(h,'EdgeColor','none')
plot(uncondPoints,krige,'color',[1 0 0],'LineWidth',2)
scatter(condPoints,condVals,'MarkerEdgeColor',[1 0 0],'LineWidth',2,'Marker','o','SizeData',50,'LineWidth',2,'MarkerFaceColor',[1,1,1]);
hold off
grid on
box on
