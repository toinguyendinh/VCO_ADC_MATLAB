%--------------Overall system----------------
Fs = 0.5*10^9;                          %sample rate
t_s = (2/Fs):(1/Fs):(101/Fs);           %cycle time
Fs_dis = 80*Fs*2^4;                     %sample rate for display
sampling_time = 80*2^4;
t = (1/Fs):(1/Fs_dis):(101/Fs);         % time range from cycle 1 to cycle 100
%---------------Input signal Block-------------------------------------------------------------%

%set up spectrum of input signal
f1 = 0.1*10^9;      %Hz
f2 = 90*10^6;       %Hz
f3 = 80*10^6;       %Hz
phi = pi*randi([0 100], 1)/50;           %rad
%Input Signal
Vpp = 0.15/7;       %scale amplitude signal
in_sig = Vpp*(3*cos(2*pi*f1*t+phi)+ cos(2*pi*f2*t+phi) + 3*cos(2*pi*f3*t+phi))+0.65;  %(V)input signal x(t) in continous time  

figure(1);    %figure out input signal and sampling frame
plot(t, in_sig);  %plot input signal x(t) in continous time
hold on;
stem(t_s, ones(1, length(t_s)));
hold off;
V_in = Vpp*(3*cos(2*pi*f1*t_s+phi)+ cos(2*pi*f2*t_s+phi) + 3*cos(2*pi*f3*t_s+phi))+0.65;
axis tight;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----VCO operation Block----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = 0; 
a=1; b=2; c=3; fc=4;    %coefficient of VCRO polynomial

f_vco =(34.33*V_in-15.15)*10^9; %(Hz) Ideal function of VCRO
%f_vco = (a*V_in^3 + b*V_in^2 + c*V_in + fc)*10^9; %approximate function of VCRO (Hz)

figure(2);  %plot ring oscillator freq according input signal Voltage
plot(V_in, f_vco); 
axis tight;
%function of VCO
VCO = [];
f_vco_con = [];
for i = 1:1:length(t_s)
    VCO = [VCO, square(2*pi*f_vco(i)*t(sampling_time*(i-1)+1:sampling_time*i))]; %Voltage signal of VCO
    f_vco_con = [f_vco_con, f_vco(i)*ones(1,sampling_time)];    % Frequency of VCO in continous time
end
VCO = 0.4*[VCO, 1];
f_vco_con = [f_vco_con, 1];

figure(3); 
plot(t, in_sig);    %plot input signal x(t) in continous time  
hold on;
plot(t, VCO);   %plot output signal of VCO 
stem(t_s, ones(1, length(t_s)));    %99 sampling frame 
hold off;
axis tight;
%--------------Intergral--------------
phase = f_vco_con(1)*(1/Fs_dis); %initial phase of VCO
for i = 2:length(f_vco_con)
    phase(i) = phase(i-1)+ f_vco_con(i)*(1/Fs_dis); %continous phase of VCO
end
figure(4)
plot(t,phase);
%-----------Counter------------------

quan_bit = 8;  %quantization bit
count = zeros(1, length(f_vco)); %set up counter
for i = 1:length(f_vco)-1
   for j = (i-1)*sampling_time+1: i*sampling_time
      if((VCO(j)<VCO(j+1)) || (j==1)) 
        count(i) = count(i)+1;
      end
      %reset counter
      if(count(i)==2^quan_bit)
          count(i)=0;
      end
   end
   count(i+1) = count(i);
end

%--------------quantizer--------------
for i = 1:1:(length(t_s)+1)
    phase_q(i) = phase(sampling_time*(i-1)+1); %quantize Phase in discrete time
end
%-----------differential block-------------%
%equivalent 

for i = 2:1:length(phase_q)
    out_sig_q(i-1) = int8(phase_q(i)-phase_q(i-1)); %output signal y(n) in discrete time
end 
plot(t_s, out_sig_q);
hold on;
grid on;
in_sig_plot = 70*in_sig - 50;
plot(t, in_sig_plot);
hold off;

%concept realization
dig_out_sig = [];
dig_out_sig(1) = count(1);  %setup differential
%decode counter
for i = 2:length(count)-1
   dig_out_sig(i) = count(i) - count(i-1);
   
   if(count(i) < count(i-1))    %decode reset counter
       dig_out_sig(i) = 2^quan_bit+count(i)-count(i-1);
   end
end

dig_out_sig_ts = [dig_out_sig, 0];
dig_out_sig_ts = dig_out_sig_ts + 20*ones(1, length(dig_out_sig_ts)) ;
figure(4);
hold on;
plot(t_s, dig_out_sig_ts);
hold off;
axis tight;

dig_out_sig = dig_out_sig - 4*ones(1, length(dig_out_sig)) ;

%-----------output to bit------------
decima_bit = out_sig_q - int8(3*ones(1, 100));
bit_stream = de2bi(decima_bit, 5, 'left-msb');
bit_stream_1 = reshape(bit_stream, [1,500]);

aaaaaa = min(dig_out_sig);

