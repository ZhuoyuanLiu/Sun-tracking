% subplot(3, 1, 1);
% plot(1:365, azi_ele_output_energy_day, '-r', 'DisplayName', 'output');  % 0th case 
% hold on;
% plot(1:365, azi_ele_output_energy_day_1st, '-r', 'DisplayName', 'output');  % 1st case
% plot(1:365, azi_ele_output_energy_day_2nd, '-r', 'DisplayName', 'output');  % 2nd case 
% xlabel('day');
% ylabel('Energy (KWH)');
% title('azi_ele Output energy');
% legend show;
% grid on;
% subplot(3, 1, 2);
% plot(1:365, azi_ele_consumed, '-b', 'DisplayName', 'Consumed');
% hold on;
% plot(1:365, azi_ele_consumed_1st, '-b', 'DisplayName', 'Consumed');
% plot(1:365, azi_ele_consumed_2nd, '-b', 'DisplayName', 'Consumed');
% xlabel('day');
% ylabel('Energy (KWH)');
% title('azi_ele Consumed energy');
% legend show;
% grid on;
% subplot(3, 1, 3);
figure,
plot(1:365, output_energy_day_0th-consumed_energy_day_0th, '-b', 'DisplayName', 'Net-0th');
hold on;
plot(1:365, output_energy_day_1st-consumed_energy_day_1st, '-r', 'DisplayName', 'Net-1st');
plot(1:365, output_energy_day_2nd-consumed_energy_day_2nd, '-g', 'DisplayName', 'Net-2nd');
plot(1:365, azi_ele_output_energy_day'-azi_ele_consumed, '-m', 'DisplayName', 'Net-4th&5th');
xlabel('Day');
ylabel('Net energy (KWH)');
title('Net energy of 4 cases');
legend show;
grid on;