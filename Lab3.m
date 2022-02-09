%% BME 517 - Lab 3
% Daniel Najarian

%% Part 2
data = csvread('lab3_part2_output.csv');

x = sqrt(data(:,1).^2+data(:,2).^2+data(:,3).^2);

figure()
plot(x,data(:,4))
xlabel('Position (mm)')
ylabel('Voltage (mV)')

disp('Intuitively, I know that the relationship sould be 1/r. That said, the graph does not reflect this, even if voltage decreases with position.')
disp('The issue lies within my COMSOL model. It failed as I could not properly delete the excess surfaces (otherwise I would have 0 current, making this trivial).')

%% Part 3
data2 = csvread('lab3_part3_output.csv');

x2 = sqrt(data2(:,1).^2+data2(:,2).^2+data2(:,3).^2);

figure()
plot(x2,data2(:,4))
xlabel('Position (mm)')
ylabel('Voltage (mV)')

new_data = sortrows(data2);
act_fxn = diff(new_data(:,4),2);

figure()
plot(x2(2:end-1),act_fxn)
xlabel('Position (mm)')
ylabel('Voltage (mV)')

disp('The activating function indicates periods of depolarization and hyperpolarization, ');
disp('so in DBS it can help researchers know if their stimulations are significant in the brain.')
