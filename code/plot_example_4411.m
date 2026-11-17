function plot_example_4411(data)
% Plot six profiles calculated at the specified steam fractions.
% The profiles are computed using one reconstructed parameter set.
% Ensemble FE distributions are plotted with plot_summary_4411.m.
colors = [0 0 1; 1 0 0];
styles = {'-','--',':'};
names = {'50%','70%','90%'};
X=linspace(0,1,numel(data(1).Ce));
Y=linspace(0,1,numel(data(1).electron_current));
% Figure A: electrons and holes
figure('Visible','off'); hold on
for k=1:3
    semilogy(X,data(k).Ce, 'Color',colors(1,:), 'LineStyle',styles{k},'LineWidth',1.5);
    semilogy(X,data(k).Ch, 'Color',colors(2,:), 'LineStyle',styles{k},'LineWidth',1.5);
end
xlabel('x/L'); ylabel('C_e, C_h (mol/m^3)'); title('4411 | Electronic carriers');
legend('e 50%','h 50%','e 70%','h 70%','e 90%','h 90%'); grid on
set(gca,'YScale','log'); % Ensure logarithmic scale with hold on.
saveas(gcf,'example_electronic_carriers.png'); close(gcf)
% Figure B: protons and vacancies
figure('Visible','off'); hold on
for k=1:3
    plot(X,data(k).Cp, 'Color',colors(1,:), 'LineStyle',styles{k},'LineWidth',1.5);
    plot(X,data(k).Cv, 'Color',colors(2,:), 'LineStyle',styles{k},'LineWidth',1.5);
end
xlabel('x/L');ylabel('C_p, C_v (mol/m^3)');title('4411 | Ionic carriers');
legend('p 50%','v 50%','p 70%','v 70%','p 90%','v 90%');grid on
saveas(gcf,'example_ionic_carriers.png'); close(gcf)
% Figure C: electron and hole currents
figure('Visible','off'); hold on
for k=1:3
    plot(Y,-data(k).electron_current, 'Color',colors(1,:), 'LineStyle',styles{k},'LineWidth',1.5);
    plot(Y,-data(k).hole_current, 'Color',colors(2,:), 'LineStyle',styles{k},'LineWidth',1.5);
end
xlabel('x/L');ylabel('i_e, i_h (A/cm^2)');title('4411 | Electronic currents');
legend('e 50%','h 50%','e 70%','h 70%','e 90%','h 90%');grid on
saveas(gcf,'example_electronic_currents.png');close(gcf)
% Figure D: proton and vacancy currents
figure('Visible','off');hold on
for k=1:3
    plot(Y,-data(k).proton_current,'Color',colors(1,:),'LineStyle',styles{k},'LineWidth',1.5);
    plot(Y,-data(k).vacancy_current,'Color',colors(2,:),'LineStyle',styles{k},'LineWidth',1.5);
end
xlabel('x/L');ylabel('i_p, i_v (A/cm^2)');title('4411 | Ionic currents');
legend('p 50%','v 50%','p 70%','v 70%','p 90%','v 90%');grid on
saveas(gcf,'example_ionic_currents.png');close(gcf)
% Figure E: electrostatic potential
figure('Visible','off');hold on
for k=1:3
    plot(X,data(k).Phi,'Color',colors(1,:),'LineStyle',styles{k},'LineWidth',1.5);
end
xlabel('x/L');ylabel('Potential (V)');title('4411 | Electrostatic potential');
legend(names);grid on
saveas(gcf,'example_potential.png');close(gcf)
% Figure F: partial pressures (left log, right linear)
figure('Visible','off');
yyaxis left;hold on
for k=1:3
    semilogy(X,data(k).pO2,'Color',[0.2 0.6 0.8],'LineStyle',styles{k},'LineWidth',1.5);
    semilogy(X,data(k).pH2,'Color',[0 0 1],'LineStyle',styles{k},'LineWidth',1.5);
end
ylabel('p_{O2}, p_{H2} (atm)');ylim([1e-30,1e10]);
yyaxis right;hold on
for k=1:3
    plot(X,data(k).pH2O,'Color',[1 0 0],'LineStyle',styles{k},'LineWidth',1.5);
end
ylabel('p_{H2O} (atm)');xlabel('x/L');title('4411 | Gas-pressure profiles');
yyaxis left; set(gca,'YScale','log'); ylim([1e-30,1e10]);
yyaxis right; set(gca,'YScale','linear');
saveas(gcf,'example_gas_pressures.png');close(gcf)
end
