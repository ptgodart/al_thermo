XRD_100c_1atm
XRD_4c_1atm
XRD_4c_68atm


figure(1)
bar(XRD_highTlowP(:,1),XRD_highTlowP(:,2),.1);
xlabel('2 Theta (Degrees)')
ylabel('Intensity (Counts)')

figure(2)
bar(XRD_lowTlowP(:,1),XRD_lowTlowP(:,2),.1);
xlabel('2 Theta (Degrees)')
ylabel('Intensity (Counts)')

figure(3)
bar(XRD_lowThighP(:,1),XRD_lowThighP(:,2),.1);
xlabel('2 Theta (Degrees)')
ylabel('Intensity (Counts)')