%parameter for XSG
fs=61.44e6;
fpga_clock=1/fs;
master_reset = zeros(200000,2);
for t = 0:199999
    master_reset(t+1,1) = t;
    master_reset(t+1,1) = master_reset(t+1,1)/fs;
end
master_reset(1,2)=1;
%master_reset(2,2)=1;
%parameter for script
L=3000;
fs=61.44e6;
f=15.36e6;
L_ifft=4096;
L_cp=288;
L_sig=3000;
%generate random bitstream 3000 sample
x=randi(10,L,1)/10;
%gateway in
x_fix=real_to_fixpt(x,1,15);
%x_fix_str=string(x_fix);
%slice 16bit to 4bit
x_rs=reshape(x_fix,L,4,4);
%16QAM
x_16qam=[];
for a=1:L
    for b=0:3
        if strcmp(x_rs(a:a,4*b+1:4*b+4),'0000')
            x_16qam(b+1,a)=1/sqrt(10)+i/sqrt(10);
        elseif strcmp(x_rs(a:a,4*b+1:4*b+4),'0001')
            x_16qam(b+1,a)=1/sqrt(10)+3*i/sqrt(10);
        
        elseif strcmp(x_rs(a:a,4*b+1:4*b+4),'0010')
            x_16qam(b+1,a)=3/sqrt(10)+i/sqrt(10);
            
        elseif strcmp(x_rs(a:a,4*b+1:4*b+4),'0011')
            x_16qam(b+1,a)=3/sqrt(10)+3*i/sqrt(10);
            
        elseif strcmp(x_rs(a:a,4*b+1:4*b+4),'0100')
            x_16qam(b+1,a)=1/sqrt(10)-i/sqrt(10);
            
        elseif strcmp(x_rs(a:a,4*b+1:4*b+4),'0101')
            x_16qam(b+1,a)=1/sqrt(10)-3*i/sqrt(10);
            
        elseif strcmp(x_rs(a:a,4*b+1:4*b+4),'0110')
            x_16qam(b+1,a)=3/sqrt(10)-i/sqrt(10);
            
        elseif strcmp(x_rs(a:a,4*b+1:4*b+4),'0111')
            x_16qam(b+1,a)=3/sqrt(10)-3*i/sqrt(10);
            
        elseif strcmp(x_rs(a:a,4*b+1:4*b+4),'1000')
            x_16qam(b+1,a)=-1/sqrt(10)+i/sqrt(10);
            
        elseif strcmp(x_rs(a:a,4*b+1:4*b+4),'1001')
            x_16qam(b+1,a)=-1/sqrt(10)+3*i/sqrt(10);
            
        elseif strcmp(x_rs(a:a,4*b+1:4*b+4),'1010')
            x_16qam(b+1,a)=-3/sqrt(10)+i/sqrt(10);
            
        elseif strcmp(x_rs(a:a,4*b+1:4*b+4),'1011')
            x_16qam(b+1,a)=-3/sqrt(10)+3*i/sqrt(10);
            
        elseif strcmp(x_rs(a:a,4*b+1:4*b+4),'1100')
            x_16qam(b+1,a)=-1/sqrt(10)-i/sqrt(10);
            
        elseif strcmp(x_rs(a:a,4*b+1:4*b+4),'1101')
            x_16qam(b+1,a)=-1/sqrt(10)-3*i/sqrt(10);
            
        elseif strcmp(x_rs(a:a,4*b+1:4*b+4),'1110')
            x_16qam(b+1,a)=-3/sqrt(10)-i/sqrt(10);
            
        else x_16qam(b+1,a)=-3/sqrt(10)-3*i/sqrt(10);
        end
    end
end
x_16qam_rs=reshape(x_16qam,4*L,1);
%padding IFFT 4096
x_pad=zeros(L_ifft*L*4/L_sig,1);
%x_pad=[]; for 3000 sample
for c=0:(L*4/L_sig-1)
    x_pad(c*L_ifft+549:c*L_ifft+549+L_sig/2-1)=x_16qam_rs(c*L_sig+1:c*L_sig+L_sig/2);
    x_pad(c*L_ifft+2050:c*L_ifft+2050+L_sig/2-1)=x_16qam_rs(c*L_sig+L_sig/2+1:c*L_sig+L_sig);
end
x_pad_rs=reshape(x_pad,L_ifft*L*4/L_sig,1);
%IFFT 4096
x_ifft=[];
for d=0:(L*4/L_sig-1)
   x_ifft(d*L_ifft+1:d*L_ifft+L_ifft)=ifft(fftshift(x_pad_rs(d*L_ifft+1:d*L_ifft+L_ifft)),L_ifft);
end
x_ifft_rs=reshape(x_ifft,(L*4/L_sig)*L_ifft,1);

%Add Cyclic Prefix
x_cp=[];
for e=0:(L*4/L_sig-1)
    x_cp(e*(L_cp+L_ifft)+1:e*(L_cp+L_ifft)+L_cp)=x_ifft_rs((e*L_ifft+L_ifft-L_cp+1):e*L_ifft+L_ifft);
    x_cp((e+1)*L_cp+e*L_ifft+1:(e+1)*L_cp+e*L_ifft+L_ifft)=x_ifft_rs(e*L_ifft+1:e*L_ifft+L_ifft);
end
x_cp_rs=reshape(x_cp,(L*4/L_sig)*(L_ifft+L_cp),1);

%plot subcarrier
subplot(5,1,1);
stem(abs(x_pad_rs(1:4096)));
xlim([1 4096]);
%plot subcarrier with shiftfft
subplot(5,1,2)
stem(abs(x_ifft(1:4096)));
xlim([1 4096]);
%plot spectrum computed ifft
subplot(5,1,3);
plot(10*log10(fftshift(abs(fft(x_ifft)))),'k');
%plot spectrum added cyclic prefix
subplot(5,1,4);
plot(10*log10(fftshift(abs(fft(x_cp)))),'k');
% plot spectrum
[PowerSpectrum,W] = pwelch(x_cp,[],[],L_ifft,fs);
subplot(5,1,5);
plot([-L_ifft/2:L_ifft/2-1]*fs/L_ifft,10*log10(fftshift(PowerSpectrum)),'k');