clc;
clear;
close all;

Jx=1;
Jy=1;
Jz=1;
Nx=100;
Ny=150;
Edge='zigzag';%edge type, 'zigzag' or 'armchair'
T=0.05;
freqoffset=0.0000000000001;
freqRange=(-0.3:0.0005:1.3)+freqoffset;%interested frequency range
VL=0.3;
VR=0;
coupling=0.1;

H=zeros(Nx*2,Nx*2,Ny);
for Nycnt=1:Ny
    if strcmpi(Edge,'armchair')
        ky=2*pi/(3*Ny)*Nycnt;
        HBA=toeplitz([2j*Jz*exp(-1j*ky) 2j*Jy*exp(1j*ky/2) zeros(1,Nx-2)],[2j*Jz*exp(-1j*ky) 2j*Jx*exp(1j*ky/2) zeros(1,Nx-2)]);
    else
        if strcmpi(Edge,'zigzag')
            ky=2*pi/(sqrt(3)*Ny)*Nycnt;
            HBA=toeplitz([2j*Jx*exp(-1j*ky*sqrt(3)/2)+2j*Jy*exp(1j*ky*sqrt(3)/2) 2j*Jz zeros(1,Nx-2)],[2j*Jx*exp(-1j*ky*sqrt(3)/2)+2j*Jy*exp(1j*ky*sqrt(3)/2)  zeros(1,Nx-1)]);
        end
    end
    HAB=HBA';
    H(:,:,Nycnt)=[zeros(Nx,Nx),HAB;HBA,zeros(Nx,Nx)];
    [V,D]=eig(H(:,:,Nycnt));
    A((Nycnt-1)*2*Nx+1:Nycnt*Nx*2)=abs(V(Nx+1,:)).^2/Ny;
    E((Nycnt-1)*2*Nx+1:Nycnt*Nx*2)=diag(D);
end

tmp1=repmat(E,length(freqRange),1); %E
tmp2=repmat(freqRange.',1,length(E))-tmp1; %omega-E

LDOS=(-pi*coupling*(tanh(tmp1/(2*T))+coth((tmp2-VL)/(2*T))).*(tmp2-VL))*A.';
RDOS=(-pi*coupling*(tanh(tmp1/(2*T))+coth((tmp2-VR)/(2*T))).*(tmp2-VR))*A.';
save('firstorderSelfE.mat','LDOS','RDOS')