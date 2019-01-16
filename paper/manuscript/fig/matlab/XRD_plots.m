XRD_100c_1atm
XRD_4c_1atm
XRD_4c_68atm
XRD_230c_68atm

figure(1)
bar(XRD_highTlowP(:,1),XRD_highTlowP(:,2),.1);
xlabel('2 Theta (Degrees)')
xlim([15 90])
ylabel('Intensity (Counts)')
title('E2: 100 C, 0.1 MPa')

figure(2)
bar(XRD_lowTlowP(:,1),XRD_lowTlowP(:,2),.1);
xlabel('2 Theta (Degrees)')
xlim([15 90]);
ylabel('Intensity (Counts)')
title('E1: 4 C, 0.1 MPa')

figure(3)
bar(XRD_lowThighP(:,1),XRD_lowThighP(:,2),.1);
xlabel('2 Theta (Degrees)')
xlim([15 90]);
ylabel('Intensity (Counts)')
title('E3: 4 C, 6.9 MPa')

figure(4)
bar(XRD_highThighP(:,1),XRD_highThighP(:,2),.1);
xlabel('2 Theta (Degrees)')
xlim([15 90]);
ylabel('Intensity (Counts)')

figure(5)
bar(XRD_highThighP_manual(:,1),XRD_highThighP_manual(:,2),.1);
xlabel('2 Theta (Degrees)')
xlim([15 90]);
ylabel('Intensity (Counts)')
title('E4: 230 C, 6.9 MPa')
