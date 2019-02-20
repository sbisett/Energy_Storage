Thickness = [0,0.03,0.04,0.05,0.08,0.1,0.12,0.15,0.18,0.2,0.3,0.4,0.5];
Loss = [13.63,3.7,3.05,2.6,1.83,1.54,1.34,1.13,0.99,0.91,0.67,0.55,0.47];
plot(Thickness,Loss)
xlabel('Insulation Thickness [m]')
ylabel('Average Percent Loss per Hour [%/hr]')
title('Insulation vs. Loss Rate in Steam Accumulator')