load ClosedLoopLPV-time.mat
load SimulinkLPV-time.mat
load linearTime.mat
figure('Name','CPU Time','NumberTitle','off');
plot(timeCL,'Linewidth',2);
hold on
plot(timeSimulink,'Linewidth',2);
plot(timeLPV_CL,'Linewidth',2);
plot(timeLPVSimulink,'Linewidth',2);
set(gca,'FontWeight','bold')
xlim([0 3000]);ylim([0 0.5])
legend('Linear (Closed-loop)','Linear (Simulink)','LPV (Closed-loop)','LPV (Simulink)','FontSize',16,'FontWeight','bold')
ylabel('Computation Time (seconds)','FontSize',16,'FontWeight','bold');xlabel('k','FontSize',16,'FontWeight','bold');
title('CPU Time','FontSize',16,'FontWeight','bold');

avgLPVSimulink = mean(timeLPVSimulink)
avgLPVCL = mean(timeLPV_CL)
avgSimulink = mean(timeSimulink)
avgCL = mean(timeCL)
