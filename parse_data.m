% load('/Users/carlosbocanegra/Downloads/data/uas4_2018-07-29-17-10-38_6.mat')

close all;

colorList = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};

%%
[Y1,X1] = ecdf(velocity_uas4(:,2));
meanX = mean(velocity_uas4(:,2));
varX = var(velocity_uas4(:,2));
[Y2,X2] = ecdf(velocity_uas4(:,3));
meanY = mean(velocity_uas4(:,3));
varY = var(velocity_uas4(:,3));
[Y3,X3] = ecdf(velocity_uas4(:,4));
meanZ = mean(velocity_uas4(:,4));
varZ = var(velocity_uas4(:,4));

figure;
subplot(1,3,1); hold on; grid minor
plot(X1,Y1,'LineWidth',2,'color',colorList{1});
[~,idx_Y1_min] = min(abs(0.2 - Y1));
[~,idx_Y1_max] = min(abs(0.8 - Y1));
plot([X1(idx_Y1_min) X1(idx_Y1_min)],ylim,'color',colorList{1},'LineStyle','--','LineWidth',1.1);
plot([X1(idx_Y1_max) X1(idx_Y1_max)],ylim,'color',colorList{1},'LineStyle','--','LineWidth',1.1);
xlabel('Velocity (m/s)','FontSize',12);
ylabel('F(x)','FontSize',12);
title('X axis','FontSize',14);

subplot(1,3,2); hold on; grid minor
plot(X2,Y2,'LineWidth',2,'color',colorList{2});
[~,idx_Y2_min] = min(abs(0.2 - Y2));
[~,idx_Y2_max] = min(abs(0.8 - Y2));
plot([X2(idx_Y2_min) X2(idx_Y2_min)],ylim,'color',colorList{2},'LineStyle','--','LineWidth',1.1);
plot([X2(idx_Y2_max) X2(idx_Y2_max)],ylim,'color',colorList{2},'LineStyle','--','LineWidth',1.1);
xlabel('Velocity (m/s)','FontSize',12);
ylabel('F(x)','FontSize',12);
title('Y axis','FontSize',14);

subplot(1,3,3); hold on; grid minor
plot(X3,Y3,'LineWidth',2,'color',colorList{3});
[~,idx_Y3_min] = min(abs(0.2 - Y3));
[~,idx_Y3_max] = min(abs(0.8 - Y3));
plot([X3(idx_Y3_min) X3(idx_Y3_min)],ylim,'color',colorList{3},'LineStyle','--','LineWidth',1.1);
plot([X3(idx_Y3_max) X3(idx_Y3_max)],ylim,'color',colorList{3},'LineStyle','--','LineWidth',1.1);
xlabel('Velocity (m/s)','FontSize',12);
ylabel('F(x)','FontSize',12);
title('Elevation','FontSize',14);

%%

[Y1,X1] = ecdf(linear_accel_uas4(:,2));
meanX = mean(linear_accel_uas4(:,2));
varX = var(linear_accel_uas4(:,2));
[Y2,X2] = ecdf(linear_accel_uas4(:,3));
meanY = mean(linear_accel_uas4(:,3));
varY = var(linear_accel_uas4(:,3));
[Y3,X3] = ecdf(linear_accel_uas4(:,4));
meanZ = mean(linear_accel_uas4(:,4));
varZ = var(linear_accel_uas4(:,4));

figure;
subplot(1,3,1); hold on; grid minor
plot(X1,Y1,'LineWidth',2,'color',colorList{1});
[~,idx_Y1_min] = min(abs(0.2 - Y1));
[~,idx_Y1_max] = min(abs(0.8 - Y1));
plot([X1(idx_Y1_min) X1(idx_Y1_min)],ylim,'color',colorList{1},'LineStyle','--','LineWidth',1.1);
plot([X1(idx_Y1_max) X1(idx_Y1_max)],ylim,'color',colorList{1},'LineStyle','--','LineWidth',1.1);
xlabel('Velocity (m/s)','FontSize',12);
ylabel('F(x)','FontSize',12);
title('X axis','FontSize',14);

subplot(1,3,2); hold on; grid minor
plot(X2,Y2,'LineWidth',2,'color',colorList{2});
[~,idx_Y2_min] = min(abs(0.2 - Y2));
[~,idx_Y2_max] = min(abs(0.8 - Y2));
plot([X2(idx_Y2_min) X2(idx_Y2_min)],ylim,'color',colorList{2},'LineStyle','--','LineWidth',1.1);
plot([X2(idx_Y2_max) X2(idx_Y2_max)],ylim,'color',colorList{2},'LineStyle','--','LineWidth',1.1);
xlabel('Velocity (m/s)','FontSize',12);
ylabel('F(x)','FontSize',12);
title('Y axis','FontSize',14);

subplot(1,3,3); hold on; grid minor
plot(X3,Y3,'LineWidth',2,'color',colorList{3});
[~,idx_Y3_min] = min(abs(0.2 - Y3));
[~,idx_Y3_max] = min(abs(0.8 - Y3));
plot([X3(idx_Y3_min) X3(idx_Y3_min)],ylim,'color',colorList{3},'LineStyle','--','LineWidth',1.1);
plot([X3(idx_Y3_max) X3(idx_Y3_max)],ylim,'color',colorList{3},'LineStyle','--','LineWidth',1.1);
xlabel('Velocity (m/s)','FontSize',12);
ylabel('F(x)','FontSize',12);
title('Elevation','FontSize',14);