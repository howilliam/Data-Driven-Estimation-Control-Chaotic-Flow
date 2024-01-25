function [PSD] = spectrum_analyser(signal,want_plot)
%spectrum_analyser
% Computes the PSD of an input signal
% -------------------------------------------------------------------------
% INPUTS:
% signal = vector (n rows, 2 columns) | Time | Magnitude |
% want_plot = (0,1) ==> (no,yes)

% OUTPUTS:
% PSD = matrix (m rows, 2 columns)

%% Signals specs
t = signal(:,1);          % Time
x = signal(:,2);          % Magnitude of signal
L = length(t);            % Length of signal  
dt = (t(end)-t(1))/(L-1); % Sampling period 
Fs = 1/dt;                % Sampling frequency   

%% FFT
X = fft(x);               % FFT
P2 = abs(X/L);            % 2-sided spectrum
P1 = P2(1:L/2+1);         % 1-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;       % Frequency range
PSD = [f; P1'];           % Output

%% Plots
if want_plot ==  1
    figure('name','Signal in time domain and PSD in frequency domain')
    
    % time domain
    subplot(2,1,1);
    plot(t,x,'k','Linewidth',1);
    title('Signal','Interpreter','Latex','Fontsize',14);
    xlabel('Time $(s)$','Interpreter','Latex','Fontsize',12);
    xlim([t(1) t(end)]);
    
    % frequency domain
    subplot(2,1,2);
    plot(PSD(1,:),PSD(2,:),'k','Linewidth',1);
    title('PSD','Interpreter','Latex','Fontsize',14);
    xlabel('$f \: (Hz)$','Interpreter','Latex','Fontsize',12);
    xlim([0 0.5]);
end

end