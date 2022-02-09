
%% BME 517 - Lab 4
% Daniel Najarian

%% Part 1

% Question 1
load('currents_big.mat')

% Question 2
r1 = [0 50 0];

% Question 3 + 4
V_ext_1 = calcVext(currents, XYZ, r1);

% Question 5
r2 = [0 100 0];
r3 = [0 200 0];
r4 = [0 300 0];

V_ext_2 = calcVext(currents, XYZ, r2);
V_ext_3 = calcVext(currents, XYZ, r3);
V_ext_4 = calcVext(currents, XYZ, r4);

t = 1:758;

subplot(2,2,1)
plot(t,V_ext_1)
title('50 um')
subplot(2,2,2)
plot(t,V_ext_2)
title('100 um')
subplot(2,2,3)
plot(t,V_ext_3)
title('200 um')
subplot(2,2,4)
plot(t,V_ext_4)
title('300 um')

% Question 6
disp('No longer detectable at 300 um distance')

%% Question 7

figure()
% Question 1
load('currents_small.mat')

V_ext_12 = calcVext(currents, XYZ, r1);
V_ext_22 = calcVext(currents, XYZ, r2);
V_ext_32 = calcVext(currents, XYZ, r3);
V_ext_42 = calcVext(currents, XYZ, r4);

subplot(2,2,1)
plot(t,V_ext_12)
title('50 um')
subplot(2,2,2)
plot(t,V_ext_22)
title('100 um')
subplot(2,2,3)
plot(t,V_ext_32)
title('200 um')
subplot(2,2,4)
plot(t,V_ext_42)
title('300 um')

disp('No longer detectable at 100 um distance')

%% Part 2

% Question 1
load('currents_big.mat')

origin = [0 0 0];
V_ext_elec = zeros(10,758);
centers = rand(10,3);
centers = centers * 50 - 25;

% Question 2
f = cell(10,1);
f{1} = 0:5^-1:1;
f{2} = 1/9:9^-1:8/9;
f{3} = 1/14:14^-1:13/14;
f{4} = 1/19:19^-1:18/19;
f{5} = 1/24:24^-1:23/24;
f{6} = 1/29:29^-1:28/29;
f{7} = 1/34:34^-1:33/34;
f{8} = 1/39:39^-1:38/39;
f{9} = 1/44:44^-1:43/44;
f{10} = 1/49:49^-1:48/49;

% Question 3
for i = 1:10
    V_ext_elec(i,:) = calcVext(currents, XYZ, centers(i,:));
end

V_ext_pt2 = mean(V_ext_elec,2);

% Question 4
tTot = 1;
dt = 0.025e-3;
cn = dsp.ColoredNoise('Color','pink','SamplesPerFrame',tTot/dt);
recordedNoise = cn();

%% Question 5
new_centers = [];
f_dead = cell(10,1);
j = 1;
for i = 1:10
    if (sqrt(centers(i,1).^2+centers(i,2).^2+centers(i,3).^2) >= 25)
        new_centers(j,:) = [centers(i,1) centers(i,2) centers(i,3)];
        f_dead{j} = f{i};
        j = j + 1;
    end
end

% Question 6
V_ext_original = cell(10,1);
figure()
subplot(1,2,1)
hold on
for i = 1:10
    arr_len = length(f{i});
    V_ext_original{i} = zeros(arr_len,1);
    for j = 1:arr_len
        V_ext_original{i}(j) = V_ext_pt2(i)+(recordedNoise(j)*1e-6);
    end
    bar(f{i},V_ext_original{i},0.1,'k')
end
title('Original')
xlabel('Time (ms)')
ylabel('Voltage (mV)')


subplot(1,2,2)
V_ext_dead = cell(size(new_centers,1),1);
for i = 1:size(new_centers,1)
    arr_len = length(f_dead{i});
    V_ext_dead{i} = zeros(arr_len,1);
    for j = 1:arr_len
        V_ext_dead{i}(j) = V_ext_pt2(i)+(recordedNoise(j)*1e-6);
    end
    bar(f_dead{i},V_ext_dead{i},0.1,'k')
end
title('Ignoring Dead Zones')
xlabel('Time (ms)')
ylabel('Voltage (mV)')

%% Part 3

% Question 1
k_matrix = csvread('lab3_part2_output.csv');
k_matrix(:,1) = k_matrix(:,1).*1e3;
k_matrix(:,2) = k_matrix(:,2).*1e3;
k_matrix(:,3) = k_matrix(:,3).*1e3;

% Question 2
load('currents_big.mat')
XYZ(:,2) = XYZ(:,2) + 50;

% Question 3
k_matrix_inter = griddata(k_matrix(:,1),k_matrix(:,2),k_matrix(:,3), ...
                     k_matrix(:,4),XYZ(:,1),XYZ(:,2),XYZ(:,3));
TF = isnan(k_matrix_inter);
k_matrix_inter(TF) = 0;

% Question 4
V_ext_electrode = currents(:,t)'*k_matrix_inter;
plot(t,V_ext_electrode);
xlabel('Time (ms)')
ylabel('Voltage (mV)')

% Question 5
disp('The shape remains similar but the magnitude has change from the order of 10^-5 to 10^-7. Note that this difference may also be due to my erraneous COMSOL model.')

%% Part 4

% Question 1
disp('33 nA is the threshold current necessary to generate action potential using a 0.1 ms duration.')

% Question 2
I_0 = ones(1,1001);
x = -2.5:0.005:2.5;
sigma = 5/(pi*10e-3^2*200e1/2000);

V_ext_neuron = I_0./(4*pi*sigma*sqrt(x.^2+1^2));

figure()
plot(x,V_ext_neuron);
xlabel('Position from end of neuron (mm)')
ylabel('Voltage (mV)')

% Question 3
g_a = 3e3;
I_intra = g_a*diff(V_ext_neuron,2);
figure()
plot(x(2:end-1),I_intra)
xlabel('Position from end of neuron (mm)')
ylabel('Voltage (mV)')

%% Question 4

I_0 = ones(1,1001);
I_0 = I_0*28e-3;
V_ext_neuron = I_0./(4*pi*sigma*sqrt(x.^2+1^2));

figure()
plot(x,V_ext_neuron);
xlabel('Position from end of neuron (mm)')
ylabel('Voltage (mV)')

% Question 3
I_intra = g_a*diff(V_ext_neuron,2);
figure()
plot(x(2:end-1),I_intra)
xlabel('Position from end of neuron (mm)')
ylabel('Voltage (mV)')


