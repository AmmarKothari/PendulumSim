%==========================================================================
% I copied most of this from something I found online.
% https://www.mathworks.com/matlabcentral/fileexchange/46101-simple-pendulum-simulation
% I am not exactly sure what we need, but I can adjust it as we figure it
% out
%==========================================================================

clear % Clears command history
clc   % Clears command window
clf  % Clears figure window
close all


PlotOn = false;
%========= Sets initial parameters for pendulum ===========================
g = 9.81;  % Gravity (ms-2) set to be upward to negative signs later on
l = 4;  % pendulum length (m)
alpha = pi/2;  % Initial angle 1
alpha_dot = 0;   % Initial angle 2
pendulummass = 1; %mass of entire pendulum
mu_alpha = 1; %friction at joint
torque_input = 1; %input torque

%====== Sets x and y coordinates of pendulum top  =========================
pendulumtopx = 0;
pendulumtopy = l;

fprintf('Single pendulum simulation for our SDM project \n')

pausetime = 0.0001;  % Pauses animation for this time
runtime = 30;  % Runs simulations for this time
tx = 0;  % Ensures time series plot remains in the figure window
time_step = 10^-1; %time step to re-evalute simulation over
steps_per_eval = 5; %number of points to evaluate over per function call


%solve equation for each time step
t_record = [];
sol_record = [];
cartesian_record = [];
torque_record = [];
t_window = linspace(0, time_step, steps_per_eval);
t_start = 0;
t_end   = time_step;
while t_end < runtime
    iterations = 1; % Sets initial iteration count to 1
    %============= Simple Controller ==================
    % i don't know if this is the best way to do this or if there is a
    % better method for inputing command torque
    K = 0.1;
%     torque_input = @(x) K.*x(2).*(pi/2 - mod(x(1), 2*pi));
    torque_input = @(x) sin(x(1));
    
    
    
    %============== Solves simple pendulum differential equations =============
    deq1=@(t,x) [x(2); -g/l * sin(x(1)) - x(2)*mu_alpha/pendulummass/l + torque_input(x)/l]; % Pendulum equations uncoupled
    [t,sol] = ode45(deq1,t_window,[alpha alpha_dot]);  % uses a numerical ode solver
%     sol = mod(sol,2*pi);
    sol_record = [sol_record; sol];
    t = t + t_start; %shift time for this time step
    t_record = [t_record; t];
    for i = 1:length(t)
        torque_record = [torque_record, torque_input(sol(i,:))];
    end
    t_start = t_end; %increment time for next loop
    t_end   = t_end + time_step;
    sol1 = sol(:,1)'; % takes the transpose for plots
    sol2 = sol(:,2)';
    alpha = sol(end,1); % update values to use in next iteration
    alpha_dot = sol(end,2);
    
    arraysize = size(t);  % Defines array size of time intervals
    timestep = t(end) - t(end-1);  % Calculates the time step of these intervals
    cartesianx = l*sin(sol1);  % Converts angles into cartesian coordinates
    cartesiany = l - l*cos(sol1);
    cartesian_record = [cartesian_record; cartesianx', cartesiany'];
    
    %============== plots results at each time interval =======================
    if PlotOn
        choice = 2; %visualize position, velocity vs. time
        for i = 1 : max(arraysize)
            subplot(2,1,1)
            plotarrayx = [pendulumtopx cartesianx(iterations)];
            plotarrayy = [pendulumtopy cartesiany(iterations)];
            plot(cartesianx(iterations),cartesiany(iterations),'ko',plotarrayx,plotarrayy,'r-')
            title(['Simple pendulum simulation            \theta = ' num2str(sol1(iterations))],'fontsize',12)
            xlabel('x [m]','fontsize',12)
            ylabel('y [m]','fontsize',12)
            f_win = 1.01;
            axis([ -l*f_win l*f_win 0 2*l*f_win])
            axis square
            subplot(2,1,2)

            % Plots either a phase portrait or time series depending on choice
            if choice == 1
                plot(sol1(iterations),sol2(iterations),'bo')
                hold on
                title('Simple pendulum phase portrait','fontsize',12)
                xlabel('\alpha','fontsize',12)
                ylabel('\alpha_dot','fontsize',12)
                axis([min(sol1) max(sol1) min(sol2) max(sol2)])
            elseif choice == 2
                plot(t(iterations),sol1(iterations),'bo')
                title(['Simple pendulum time series for \alpha1       t = ' num2str(t(iterations))],'fontsize',12)
                xlabel('t [seconds]','fontsize',12)
                ylabel('\alpha1','fontsize',12)
                hold on  % Holds previous values
                axis([0 t(iterations)+(t(2)-t(1)) min(sol_record(:,1)) max(sol_record(:,1))])
                tx = tx + timestep;  % Aligns results with the figure window
            end
            pause(pausetime)  % Shows results at each time interval
            iterations = iterations + 1;  % increases iteration count by 1
        end
    end
    
    
end
% l = sqrt(sum((cartesian_record - [0,4]).^2,2));
% figure(); plot(l);

% plot(t, sol1, 'b', t, sol2, 'r');
% figure()
% plot(t, l*sin(sol1), 'b', t, l*cos(sol1), 'r');
figure()
subplot(3,1,1)
plot(t_record, sol_record(:,1))
subplot(3,1,2)
plot(cartesian_record(:,1), cartesian_record(:,2),'o')
subplot(3,1,3)
plot(t_record, torque_record)

%=========================== End of program ===============================
