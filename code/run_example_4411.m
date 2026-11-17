% Forward electrochemical transport example for BZCYYb4411 at 600 C.
% Requires MATLAB Optimization Toolbox (fsolve).
% Uses one reconstructed material-specific parameter set and -1 A/cm2.
base_dir = fileparts(mfilename('fullpath'));
old_dir = pwd;
cd(base_dir);
try
    fractions = [0.50, 0.70, 0.90]; % Oxygen-electrode steam fractions.
    hydrogen_fraction = 0.99; % H2 mole fraction at the hydrogen electrode.
    oxygen = 0.2*(1-fractions);  % 0.10 / 0.06 / 0.02 atm
    for k=1:numel(fractions)
        [yinit,CeS,ChS,CpS,CvS,CoS,PhiS,CoSOCV,PO2SOCV,PO2S,PH2S,PH2OS, ...
         rleak,rleakS,rv,rp,re,rh,icS,V] = solve_transport_steady( ...
          -1,600,1e-5,0.05,0.05,0.005,hydrogen_fraction,fractions(k),oxygen(k),0.19,0,1);
        data(k).fraction=fractions(k);
        data(k).Ce=CeS;
        data(k).Ch=ChS;
        data(k).Cp=CpS;
        data(k).Cv=CvS;
        data(k).Phi=PhiS;
        data(k).pO2=PO2S;
        data(k).pH2=PH2S;
        data(k).pH2O=PH2OS;
        data(k).electron_current=re;
        data(k).hole_current=rh;
        data(k).proton_current=rp;
        data(k).vacancy_current=rv;
        data(k).FE=rp(100);
        fprintf('4411 | %.0f%% H2O | FE = %.6f%% at -1 A/cm2\n',100*fractions(k),100*rp(100));
    end
    save('example_results_4411.mat','data');
    plot_example_4411(data); % Plot the spatial concentration, current, and potential profiles.
catch ME
    cd(old_dir);
    rethrow(ME);
end
cd(old_dir);
