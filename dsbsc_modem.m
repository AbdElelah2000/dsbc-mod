clear
hold off
format long e
% DSB-SC
% Abd Elelah Arafah - arafaha - 400197623

N = 4096; %No. of FFT samples
sampling_rate = 100.0e3; %unit Hz
tstep = 1/sampling_rate;
tmax = N*tstep/2;

tmin = -tmax;
tt = tmin:tstep:tmax-tstep;
fmax = sampling_rate/2; 
fmin = -fmax;
fstep = (fmax-fmin)/N;
freq = fmin:fstep:fmax-fstep;

fc=10e3;
Ac = 1;
Am = 1;
fm = 1e3;
ct=Ac*cos(2*pi*fc*tt);
mt = Am*cos(2*pi*fm*tt);
st = mt.*ct;

%Carrier in Time domain
figure(1)
Hp1 = plot(tt,ct);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('Carrier c(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Carrier : Time domain');
axis([-1e-3 1e-3 -1.1 1.1])
pause(1)

%message signal in Time domain
figure(2)
Hp1 = plot(tt,mt);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('message  m(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('message signal : Time domain');
axis([-0.01 0.01 0 1.1])
pause(1)

%DSB-SC modulated wave in Time domain
figure(3)
Hp1 = plot(tt,st);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('s(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('DSB-SC modulated wave : Time domain');
axis([-10/fc 10/fc min(st) max(st)])
pause(1)

Mf = fftshift(fft(fftshift(mt)))/(2*fmax);
%Spectrum of the message signal in frequency domain
figure(4)
%The amplitude of the spectrum is different from the Fourier transform
%amplitude due to discretization of discrete Fourier transform
Hp1=plot(freq,abs(Mf));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|M(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum of the message signal');
axis ([-15e3 15e3 0 max(abs(Mf))])
pause(1)


Sf = fftshift(fft(fftshift(st)))/(2*fmax);

%Spectrum of the DSB-SC wave in frequency domain
figure(5)
Hp1=plot(freq,abs(Sf));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|S(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum of the DSB-SC wave S(f)');
axis ([-30e3 30e3 0 max(abs(Sf))])
pause(1)

%DSB-SC demodulation

%Local oscillator at the receiver perfectly synchronized
thet=0;
lo = cos(2*pi*fc*tt + thet); 
st1 = st .* lo;

%signal after remodulation at Rx in time domain
figure(6)
Hp1=plot(tt,st1);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('  s hat(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('signal after remodulation at Rx, s hat(t) : Time domain');
axis([-10/fc 10/fc min(st1) max(st1)])
pause(1)
Sf1 = fftshift(fft(fftshift(st1)))/(2*fmax);

% Spectrum S hat(f) in freq domain
figure(7)
Hp1=plot(freq,abs(Sf1));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|S hat(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum S hat(f)');
axis ([-50e3 50e3 0 max(abs(Sf1))])
pause(1)

%Low pass filtering
f_cutoff = 1.1*fm; %as chosen in lab 2!

%ideal low pass filter code
n=1;
for f = freq
    if abs(f) < f_cutoff
        Hf(n) = 1;
    else
        Hf(n) = 0;
    end
n=n+1;
end

%Output of low pass filter in time domain
Mf1 = Sf1 .* Hf;
mt1 = 2*fmax*fftshift(ifft(fftshift(Mf1)));
figure(8)
Hp1=plot(tt,mt1,'r',tt,mt*0.5,'g.');
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('m hat(t)(Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Output of low pass filter,  m hat(t)  : Time domain');
axis([-0.01 0.01 min(mt*0.5) max(mt*0.5)])
legend('LPF output', 'message sig');

% Abd Elelah Arafah - arafaha - 400197623













     
