clc;
clear;close all;
warning('off','all')
global T
 
Nx=70;%length along x direction
Edge='armchair';%edge type, 'zigzag' or 'armchair'
VL=0.7;%range of left spin bias
T=0.005;%temperature
Jx=1;Jy=1;Jz=1; %equals J_{ij} in Lee's paper, not J_\lambda
freqoffset=0.0000000000001;
freqstep=0.00001;
freqRange=(-0.3:freqstep:1.3)+freqoffset;%interested frequency range
kyRange=2*pi/(3*100)*(0:1:20);%range of ky
RealSelfE=0;%1 to turn on the real part of self energy, 0 to off
couplingratio=10; %coupling=couplingratio times 0.1
VR=0;%right spin bias
func=@(mu,Temp,x) 2*tanh((x-mu)/(2*Temp));

input=load('firstorderSelfE.mat');%load self energy calculated by YLselfE.m
inputrange=-0.3:0.0005:1.3;%specify the frequency range of loaded self energy

GfR=zeros(Nx*2,Nx*2);
sigfR=zeros(Nx*2,Nx*2);
trans=zeros(size(freqRange));

LsigfR=interp1(inputrange,1j*couplingratio*input.LDOS,freqRange-freqoffset);%interpolate the self energy
RsigfR=interp1(inputrange,1j*couplingratio*input.RDOS,freqRange-freqoffset);%interpolate the self energy
LsigfK=imag(LsigfR).*2.*tanh((freqRange-VLrange)/(2*T))*1j;
RsigfK=imag(RsigfR).*2.*tanh((freqRange-VR)/(2*T))*1j;
Ifreq=zeros(length(freqRange),length(kyRange));
kycnt=1;

for ky=kyRange
    if strcmpi(Edge,'armchair')
        HBA=toeplitz([2j*Jz*exp(-1j*ky) 2j*Jy*exp(1j*ky/2) zeros(1,Nx-2)],[2j*Jz*exp(-1j*ky) 2j*Jx*exp(1j*ky/2) zeros(1,Nx-2)]);
    else
        if strcmpi(Edge,'zigzag')
            HBA=toeplitz([2j*Jx*exp(-1j*ky*sqrt(3)/2)+2j*Jy*exp(1j*ky*sqrt(3)/2) 2j*Jz zeros(1,Nx-2)],[2j*Jx*exp(-1j*ky*sqrt(3)/2)+2j*Jy*exp(1j*ky*sqrt(3)/2)  zeros(1,Nx-1)]);
        end
    end
    HAB=HBA';
    H=[zeros(Nx,Nx),HAB;HBA,zeros(Nx,Nx)];
    freqcnt=1;
    for omega=freqRange
        if strcmpi(Edge,'armchair')
            sigfR(Nx,Nx)=RsigfR(freqcnt);
            sigfR(2*Nx,2*Nx)=RsigfR(freqcnt);
            sigfR(Nx+1,Nx+1)=LsigfR(freqcnt);
            sigfR(1,1)=LsigfR(freqcnt);
        else
            if strcmpi(Edge,'zigzag')
                sigfR(Nx,Nx)=RsigfR(freqcnt);
                sigfR(Nx+1,Nx+1)=LsigfR(freqcnt);
            end
        end
        
        if RealSelfE==1
            GfR=inv(eye(size(H,1))*omega-H-sigfR);
        else
            if RealSelfE==0
                GfR=inv(eye(size(H,1))*omega-H-1j*imag(sigfR));
            end
        end
        if strcmpi(Edge,'armchair')
            GfRsq=abs(GfR([1,Nx,Nx+1,Nx*2],[1,Nx,Nx+1,Nx*2])).^2;
            trans(freqcnt)=-(GfRsq(1,2)+GfRsq(1,4)+GfRsq(3,2)+GfRsq(3,4))*LsigfR(freqcnt)*RsigfR(freqcnt)*4;%transmission coefficient
        else
            if strcmpi(Edge,'zigzag')
                GfRsq=abs(GfR([1,Nx,Nx+1,Nx*2],[1,Nx,Nx+1,Nx*2])).^2;
                trans(freqcnt)=-(GfRsq(3,2))*LsigfR(freqcnt)*RsigfR(freqcnt)*4;%transmission coefficient
            end
        end
        freqcnt=freqcnt+1;
    end
    fLR=fermiD(freqRange-VLrange)-fermiD(freqRange);%fL-fR
    Ifreq(:,kycnt)=fLR.*trans;
    I(kycnt)=trapz(trans.*fLR)*freqstep/(2*pi);%current in ky space
    kycnt=kycnt+1;
end

%% Functions
function n=fermiD(E)
global T
n=1./(exp(E/T)+1);
end