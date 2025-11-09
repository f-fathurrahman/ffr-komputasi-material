%copyright by J. E Hasbun and T. Datta
% ch10_ising.m - by William D. Baez, CurtisLee Thornton, and T. Datta
% GUI implementation by -  CurtisLee Thornton and Alexander Price

% Monte Carlo GUI code is available from the publisher separately

% Ising model code

clearvars;
close all

%{
VARIABLE DEFINITIONS:
    latsize-----------(scalar) Length of lattice edges
    N-----------------(scalar) Number of sites in the Lattice
    initT-------------(scalar) Initial Temperature
    finT--------------(scalar) Final Temperature
    Tdt---------------(scalar) Temperature Step Size
    Temp--------------(array)  Temperature Array
    MCSS--------------(scalar) Monte Carlo Steps per Site
    MCS---------------(scalar) Monte Carlo Steps
    M-----------------(array)  Initial magnetization (all up)
    magnetization-----(array)  Magnetization for each MCSS
    Emagnetization----(array)  Magnetization for Equilibrium
    average_mag-------(array)  Average Magnetization
    susceptibility----(array)  Susceptibility Array
    Ediff-------------(scalar) Total Energy Contribution from Nearest-neighbor
    expect_mag--------(scalar) Expectation Value for Magnetization
    expect_magsqr-----(scalar) Expectation Value for Magnetization Squared
    mag_sumsqr--------(scalar) Sum of the Squared Magnetization
    initH-------------(scalar) Magnitude of External Magnetic Field

NOTE: Comment or uncomment as need for magnetzation and susceptibility or
equilibration study.
 %}

% % % % User Input Variables % % % %
latsize = input('Enter 2-D Lattice Dimension:');
MCSS    = input('Enter the number of Monte Carlo Steps per Site:');
initH   = input('Enter the Magnetic Field Strength:');
initT   = input('Enter the initial Temperature:');
finT    = input('Enter the final Temperature:');
Tdt     = input('Enter the Temperature Step Size:');

% % % % Generates initial lattice with all spins up % % % %
M = ones(latsize);

% % Random Initial State (USED FOR EQUILIBRATION STUDY ONLY) % % % %
% Comment this block if not running equilibration study else uncomment.
for ix = 1:latsize
    for iy = 1:latsize
        if rand < 0.5
            M(ix,iy) = 1;
        else
            M(ix,iy) = -1;
        end
    end
end

% % % % Temperature Array % % % %
Temp = initT:Tdt:finT;

% % % % Total Number of Lattice Sites % % % %
N = latsize*latsize;

% % % % Total Number of Monte Carlo Steps % % % %
MCS = MCSS*N;

% % % % Initilize Arrays % % % %
magnetization = zeros(1,MCSS);
Emagnetization = zeros(1,MCS);
average_mag = zeros(1,length(Temp));
susceptibility = zeros(1,length(Temp));

% % % % Start Stop Watch % % % %
tstart = tic;

% % % % Temperature Step % % % %
for t = 1:length(Temp)

    % % % % Prints the Current Temperature at Every Step % % % %
    fprintf('Current Temp: %g\n',Temp(t))

    % % % % Starts the counter % % % %
    count = 0;
    countE = 0;

    % % % % Montecarlo Step % % % %
    for montecarlosteps = 1:MCS

        % % % % Select a Random Lattice Site % % % %
        a = floor(rand*latsize + 1); % row number
        b = floor(rand*latsize + 1); % column number

        % % % % Energy for the Spin ABOVE the Lattice Point % % % %
        if a == 1
            Edifftop = M(latsize,b);
        else
            Edifftop = M(a-1,b);
        end

        % % % % Energy for the Spin BELOW the Lattice Point % % % %
        if a == latsize
            Ediffbottom = M(1,b);
        else
            Ediffbottom = M(a+1,b);
        end

        % % % % Energy for the Spin LEFT of the Lattice Point % % % %
        if b == 1
            Ediffleft = M(a,latsize);
        else
            Ediffleft = M(a,b-1);
        end

        % % % % Energy for the Spin RIGHT of the Lattice Point % % % %
        if b == latsize
            Ediffright = M(a,1);
        else
            Ediffright = M(a,b+1);
        end

        % % % % Total Energy Difference % % % %
        Ediff = 2*M(a,b)*(Edifftop + Ediffbottom + Ediffleft + Ediffright)...
                + (2*M(a,b)*initH);

        % % % % Spin Flip Decision % % % %
        if Ediff <= 0
            M(a,b) = -M(a,b);
        elseif rand < exp(-Ediff/Temp(t))
            M(a,b) = -M(a,b);
        end

        % % % % Magnetization Tracker per FLIP % % % %
        countE = countE + 1;
        Emagnetization(countE) = abs((sum(sum(M))))/N;

        % % % % Magnetization Tracker per SWEEP % % % %
        if mod(montecarlosteps,N) == 0
            count = count + 1;
            magnetization(count) = abs((sum(sum(M))));
        end
    end

    % % % % Calculations % % % %
    average_mag(t)= abs(((sum(magnetization(((0.90)*MCSS+1):MCSS))/(0.10*MCSS))))/N;

    expect_mag = sum(magnetization)/length(magnetization);
    expect_magsqr = expect_mag^2;

    mag_sumsqr = sum(magnetization.^2)/length(magnetization);

    susceptibility(t) = ((mag_sumsqr - expect_magsqr)/Temp(t))/N;
end

% % % % End Stop Watch % % % %
telapsed = toc(tstart);
fprintf('Total Calculation Time: %g seconds\n',telapsed)

% % % % Plots % % % %

% USED FOR MAGNETIZATION AND SUSCEPTIBILITY STUDY ONLY
% Comment block if not running equilibration study, uncommented for magnetization
% and susceptibility study
% Select the block of code below and use Ctrl + R to comment
% Select the block of code below and use Ctrl + T to uncomment

figure('Name','Plots','NumberTitle','off')
subplot 121
title('Magnetic Phase Diagram')
plot(Temp(2:length(Temp)),average_mag(2:length(Temp)),'-bo')
xlabel('Temperature (K)');
ylabel('Average Magnetization per Spin (M/N)');
subplot 122
title('Susceptibility');
plot(Temp(2:length(Temp)),susceptibility(2:length(Temp)),'-bo')
xlabel('Temperature (K)');
ylabel('Susceptibility per Spin (S/N)');

% USED FOR EQUILIBRATION STUDY ONLY
% Comment block if not running equilibration study, else uncomment
% Select the block of code below and use Ctrl + R to comment
% Select the block of code below and use Ctrl + T to uncomment

% figure('Name','Equilibration','NumberTitle','off')
% plot(Emagnetization,'-b*')
% title('Equilibration')
% ylim([0 1])
% xlabel('MCS (Monte Carlo Steps)')
% ylabel('Magnetization per Spin (M/S)')
