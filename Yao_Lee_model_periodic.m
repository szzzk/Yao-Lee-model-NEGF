clc;
clear;close all;
warning('off','all')

global Nx Ny delta T GammaL GammaR cutoff IntfreqRange freqRange freqstep
Nx=11;
Ny=9; %Number of unit cells
Jx=1;Jy=1;Jz=1; %equals J_{ij} in Lee's paper, not J_\lambda

VLrange=0.1;
freqend=7;
freqoffset=0.0000000000001;
freqstep=0.002;
freqRange=(-freqend:freqstep:freqend)+freqoffset;%The frequency range that we are interested in
IntfreqRange=(-freqend*2:freqstep:freqend*2)+freqoffset; %The freq range needed for integral, which should be at least 2 times wider than freqRange
Lcoupling=0.1;
Rcoupling=0.1;
Edge='zigzag';%edge type, 'zigzag' or 'armchair'
RealSelfE=0;%1 to turn on the real part of self energy, 0 to off
delta=0.01; %The infinitesimal appearing in the Green's functions
VR=0;
T=0.1;
cutoff=100;%cutoff set in the spectral density of the effective bath, just set to a large number

func=@(mu,Temp,x) 2*tanh((x-mu)/(2*Temp));
%% Hamiltonian construction
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
end

%% Bare Green's functions
gfR=zeros(2*Nx,2*Nx,length(freqRange),Ny);
gfA=zeros(2*Nx,2*Nx,length(freqRange),Ny);
gfK=zeros(2*Nx,2*Nx,length(freqRange),Ny);
for Nycnt=1:Ny
    if strcmpi(Edge,'armchair')
        ky=2*pi/(3*Ny)*Nycnt;
    else
        if strcmpi(Edge,'zigzag')
            ky=2*pi/(sqrt(3)*Ny)*Nycnt;
        end
    end
    [V,D]=eig(H(:,:,Nycnt));
    Energy=diag(D);
    freqcnt=1;
    for omega=freqRange
        gfR(:,:,freqcnt,Nycnt)=V*diag(1./(omega-Energy+1j*delta))*V';
        gfA(:,:,freqcnt,Nycnt)=gfR(:,:,freqcnt,Nycnt)';
        gfK(:,:,freqcnt,Nycnt)=V*diag(-2j*delta./((omega-Energy).^2+delta^2).*(1-2*fermiD(Energy)))*V';
        freqcnt=freqcnt+1;
    end
end
gfspec=(gfR-gfA)/(-2j*pi);

%% Full Green's functions

biascnt=1;
maxiter=40;

DiagGL=zeros(Nx*2,1);
DiagGR=zeros(Nx*2,1);
if strcmpi(Edge,'armchair')
    DiagGL([1 Nx+1])=Lcoupling*cutoff;
    DiagGR([Nx 2*Nx])=Rcoupling*cutoff;
else
    if strcmpi(Edge,'zigzag')
        DiagGL(Nx+1)=Lcoupling*cutoff;
        DiagGR(Nx)=Rcoupling*cutoff;
    end
end
GammaL=diag(DiagGL);
GammaR=diag(DiagGR);

midpt=(length(IntfreqRange)+1)/2;

global LbathR RbathR LbathA RbathA LbathK RbathK;
global LbathR2 RbathR2 LbathA2 RbathA2 LbathK2 RbathK2
LbathR=zeros(length(freqRange),length(freqRange));
RbathR=zeros(length(freqRange),length(freqRange));
LbathA=zeros(length(freqRange),length(freqRange));
RbathA=zeros(length(freqRange),length(freqRange));
LbathK=zeros(length(freqRange),length(freqRange));
RbathK=zeros(length(freqRange),length(freqRange));
LbathR2=zeros(length(freqRange),length(freqRange));
RbathR2=zeros(length(freqRange),length(freqRange));
LbathA2=zeros(length(freqRange),length(freqRange));
RbathA2=zeros(length(freqRange),length(freqRange));
LbathK2=zeros(length(freqRange),length(freqRange));
RbathK2=zeros(length(freqRange),length(freqRange));
GfR=zeros(Nx*2,Nx*2,length(freqRange),Ny);
GfA=zeros(Nx*2,Nx*2,length(freqRange),Ny);
GfK=zeros(Nx*2,Nx*2,length(freqRange),Ny);
GnuR=zeros(Nx*2,Nx*2,length(freqRange),Ny);
GnuA=zeros(Nx*2,Nx*2,length(freqRange),Ny);
GnuK=zeros(Nx*2,Nx*2,length(freqRange),Ny);

for VL=VLrange
    for ii=1:length(freqRange)
        LbathR(:,ii)=-pi*(VL+1j*cutoff)./(IntfreqRange(midpt+ii-1:-1:ii)+1j*cutoff)*Lcoupling*cutoff;
        RbathR(:,ii)=-pi*(VR+1j*cutoff)./(IntfreqRange(midpt+ii-1:-1:ii)+1j*cutoff)*Rcoupling*cutoff;
        LbathA(:,ii)=-pi*(VL-1j*cutoff)./(IntfreqRange(midpt+ii-1:-1:ii)-1j*cutoff)*Lcoupling*cutoff;
        RbathA(:,ii)=-pi*(VR-1j*cutoff)./(IntfreqRange(midpt+ii-1:-1:ii)-1j*cutoff)*Rcoupling*cutoff;
        LbathK(:,ii)=-2j*pi*(IntfreqRange(midpt+ii-1:-1:ii)-VL).*cutoff./(IntfreqRange(midpt+ii-1:-1:ii).^2+cutoff^2).*coth((IntfreqRange(midpt+ii-1:-1:ii)-VL)/(2*T))*Lcoupling*cutoff;
        RbathK(:,ii)=-2j*pi*(IntfreqRange(midpt+ii-1:-1:ii)-VR).*cutoff./(IntfreqRange(midpt+ii-1:-1:ii).^2+cutoff^2).*coth((IntfreqRange(midpt+ii-1:-1:ii)-VR)/(2*T))*Rcoupling*cutoff;
        
        LbathR2(:,ii)=-pi*(VL+1j*cutoff)./(IntfreqRange(midpt-ii+1:end-ii+1)+1j*cutoff)*Lcoupling*cutoff;
        RbathR2(:,ii)=-pi*(VR+1j*cutoff)./(IntfreqRange(midpt-ii+1:end-ii+1)+1j*cutoff)*Rcoupling*cutoff;
        LbathA2(:,ii)=-pi*(VL-1j*cutoff)./(IntfreqRange(midpt-ii+1:end-ii+1)-1j*cutoff)*Lcoupling*cutoff;
        RbathA2(:,ii)=-pi*(VR-1j*cutoff)./(IntfreqRange(midpt-ii+1:end-ii+1)-1j*cutoff)*Rcoupling*cutoff;
        LbathK2(:,ii)=-2j*pi*(IntfreqRange(midpt-ii+1:end-ii+1)-VL).*cutoff./(IntfreqRange(midpt-ii+1:end-ii+1).^2+cutoff^2).*coth((IntfreqRange(midpt-ii+1:end-ii+1)-VL)/(2*T))*Lcoupling*cutoff;
        RbathK2(:,ii)=-2j*pi*(IntfreqRange(midpt-ii+1:end-ii+1)-VR).*cutoff./(IntfreqRange(midpt-ii+1:end-ii+1).^2+cutoff^2).*coth((IntfreqRange(midpt-ii+1:end-ii+1)-VR)/(2*T))*Rcoupling*cutoff;

    end    
    [signuR,~,signuK]=nuselfE(sum(gfR,4)/Ny,sum(gfA,4)/Ny,sum(gfK,4)/Ny,VL,VR);
    
    for iter=1:maxiter        
        for Nycnt=1:Ny
            freqcnt=1;
            for omega=freqRange
                if RealSelfE==1
                    GnuR(:,:,freqcnt,Nycnt)=inv((eye(size(D))*omega-H(:,:,Nycnt))/(2j)-signuR(:,:,freqcnt));
                else
                    if RealSelfE==0
                        GnuR(:,:,freqcnt,Nycnt)=inv((eye(size(D))*omega-H(:,:,Nycnt))/(2j)-real(signuR(:,:,freqcnt)));
                    end
                end
                GnuA(:,:,freqcnt,Nycnt)=-GnuR(:,:,freqcnt,Nycnt)';
                GnuK(:,:,freqcnt,Nycnt)=GnuR(:,:,freqcnt,Nycnt)*signuK(:,:,freqcnt)*GnuA(:,:,freqcnt,Nycnt);
                freqcnt=freqcnt+1;
            end
        end
        [sigfR,sigfA,sigfK]=fselfE(sum(GnuR,4)/Ny,sum(GnuA,4)/Ny,sum(GnuK,4)/Ny,VL,VR);

        formerGfR=GfR;
        for Nycnt=1:Ny
            freqcnt=1;
            for omega=freqRange
                if RealSelfE==1
                    GfR(:,:,freqcnt,Nycnt)=inv(eye(size(D))*omega-H(:,:,Nycnt)-sigfR(:,:,freqcnt));
                else
                    if RealSelfE==0
                        GfR(:,:,freqcnt,Nycnt)=inv(eye(size(D))*omega-H(:,:,Nycnt)-1j*imag(sigfR(:,:,freqcnt)));
                    end
                end                
                GfA(:,:,freqcnt,Nycnt)=GfR(:,:,freqcnt,Nycnt)';
                GfK(:,:,freqcnt,Nycnt)=GfR(:,:,freqcnt,Nycnt)*sigfK(:,:,freqcnt)*GfA(:,:,freqcnt,Nycnt);
                freqcnt=freqcnt+1;
            end
        end
        [signuR,~,signuK]=nuselfE(sum(GfR,4)/Ny,sum(GfA,4)/Ny,sum(GfK,4)/Ny,VL,VR);
        iterdiff(iter)=max(max(max(max(abs(GfR-formerGfR)))));
        if iterdiff(iter)<=10^(-8)
            break
        end        
    end
    if iter==maxiter
        disp(['V=',num2str(VL),',havent converged, iterdiff=',num2str(iterdiff(maxiter))])
        %I(biascnt)=nan;
    else
        disp(['V=',num2str(VL),',have converged after itering ',num2str(iter),' times'])
    end
    Gfl=1/2*(GfK-GfR+GfA);    
    Gfspec=(GfR-GfA)/(-2j*pi);
    Sumspec=reshape(sum(Gfspec,3)*freqstep,Nx*2,Nx*2,Ny);    
    I(:,:,biascnt)=(-sum(H.*permute(reshape(sum(Gfl,3),size(H)),[2 1 3]),3)+sum(permute(H,[2 1 3]).*reshape(sum(Gfl,3),size(H)),3))*freqstep/(2*pi);        
    disfL(:,biascnt)=imag(squeeze(sigfK(Nx+1,Nx+1,:)))./imag(squeeze(sigfR(Nx+1,Nx+1,:)));%effective bath fermi distribution
    disfR(:,biascnt)=imag(squeeze(sigfK(Nx,Nx,:)))./imag(squeeze(sigfR(Nx,Nx,:)));
    disnuL(:,biascnt)=real(squeeze(signuK(Nx+1,Nx+1,:)))./real(squeeze(signuR(Nx+1,Nx+1,:)));
    disnuR(:,biascnt)=real(squeeze(signuK(Nx,Nx,:)))./real(squeeze(signuR(Nx,Nx,:)));
    coefffL(:,biascnt)=coeffvalues(fit(reshape(freqRange,[],1),disfL(:,biascnt),func));
    coefffR(:,biascnt)=coeffvalues(fit(reshape(freqRange,[],1),disfR(:,biascnt),func));
    coeffnuL(:,biascnt)=coeffvalues(fit(reshape(freqRange,[],1),disnuL(:,biascnt),func));
    coeffnuR(:,biascnt)=coeffvalues(fit(reshape(freqRange,[],1),disnuR(:,biascnt),func));
    biascnt=biascnt+1;
end


SpecErr=abs(Sumspec-repmat(eye(Nx*2,Nx*2),1,1,Ny)); %Error of spectral function integrated over frequency domain
for ii=1:size(Sumspec,3)
    offd(:,:,ii)=SpecErr(:,:,ii)-diag(diag(SpecErr(:,:,ii))); %off diagonal error
    dia(:,:,ii)=diag(SpecErr(:,:,ii));%diagonal error
end
offErr=max(max(max(offd)));
diagErr=max(max(max(dia)));


%% Functions
function [sigfR,sigfA,sigfK]=fselfE(GnuR,GnuA,GnuK)
global  GammaL GammaR Nx freqRange freqstep
global LbathR RbathR LbathA RbathA LbathK RbathK
LInd=reshape(find(diag(GammaL)~=0),1,[]);
RInd=reshape(find(diag(GammaR)~=0),1,[]);

sigfR=zeros(Nx*2,Nx*2,length(freqRange));
sigfK=zeros(Nx*2,Nx*2,length(freqRange));
sigfA=zeros(Nx*2,Nx*2,length(freqRange));
LGnuR=zeros(length(LInd),length(freqRange));
LGnuA=zeros(length(LInd),length(freqRange));
LGnuK=zeros(length(LInd),length(freqRange));
RGnuR=zeros(length(RInd),length(freqRange));
RGnuA=zeros(length(RInd),length(freqRange));
RGnuK=zeros(length(RInd),length(freqRange));

for ii=1:length(LInd)
LGnuK(ii,:)=GnuK(LInd(ii),LInd(ii),:);
LGnuR(ii,:)=GnuR(LInd(ii),LInd(ii),:);
LGnuA(ii,:)=GnuA(LInd(ii),LInd(ii),:);
end
for ii=1:length(RInd)
RGnuK(ii,:)=GnuK(RInd(ii),RInd(ii),:);
RGnuR(ii,:)=GnuR(RInd(ii),RInd(ii),:);
RGnuA(ii,:)=GnuA(RInd(ii),RInd(ii),:);
end

LsigfR=1/(4*pi)*(LGnuK*LbathR+LGnuR*LbathK)*freqstep;
RsigfR=1/(4*pi)*(RGnuK*RbathR+RGnuR*RbathK)*freqstep;
LsigfK=1/(4*pi)*1j*imag((LGnuR-LGnuA)*(LbathR-LbathA)+LGnuK*LbathK)*freqstep;
RsigfK=1/(4*pi)*1j*imag((RGnuR-RGnuA)*(RbathR-RbathA)+RGnuK*RbathK)*freqstep;
LsigfA=conj(LsigfR);
RsigfA=conj(RsigfR);

for ii=1:length(LInd)
    sigfR(LInd(ii),LInd(ii),:)=LsigfR(ii,:);
    sigfA(LInd(ii),LInd(ii),:)=LsigfA(ii,:);
    sigfK(LInd(ii),LInd(ii),:)=LsigfK(ii,:);
end
for ii=1:length(RInd)
    sigfR(RInd(ii),RInd(ii),:)=RsigfR(ii,:);
    sigfA(RInd(ii),RInd(ii),:)=RsigfA(ii,:);
    sigfK(RInd(ii),RInd(ii),:)=RsigfK(ii,:);
end
end

function [signuR,signuA,signuK]=nuselfE(GfR,GfA,GfK)
global  GammaL GammaR Nx freqRange freqstep
global LbathR2 RbathR2 LbathA2 RbathA2 LbathK2 RbathK2
LInd=reshape(find(diag(GammaL)~=0),1,[]);
RInd=reshape(find(diag(GammaR)~=0),1,[]);
signuR=zeros(Nx*2,Nx*2,length(freqRange));
signuK=zeros(Nx*2,Nx*2,length(freqRange));
signuA=zeros(Nx*2,Nx*2,length(freqRange));
LGfR=zeros(length(LInd),length(freqRange));
LGfA=zeros(length(LInd),length(freqRange));
LGfK=zeros(length(LInd),length(freqRange));
RGfR=zeros(length(RInd),length(freqRange));
RGfA=zeros(length(RInd),length(freqRange));
RGfK=zeros(length(RInd),length(freqRange));

for ii=1:length(LInd)
LGfK(ii,:)=GfK(LInd(ii),LInd(ii),:);
LGfR(ii,:)=GfR(LInd(ii),LInd(ii),:);
LGfA(ii,:)=GfA(LInd(ii),LInd(ii),:);
end
for ii=1:length(RInd)
RGfK(ii,:)=GfK(RInd(ii),RInd(ii),:);
RGfR(ii,:)=GfR(RInd(ii),RInd(ii),:);
RGfA(ii,:)=GfA(RInd(ii),RInd(ii),:);
end

LsignuR1=1/(4*pi)*(LGfK*LbathA2+LGfR*LbathK2)*freqstep;
RsignuR1=1/(4*pi)*(RGfK*RbathA2+RGfR*RbathK2)*freqstep;
LsignuK1=1/(4*pi)*real((LGfR-LGfA)*(LbathA2-LbathR2)+LGfK*LbathK2)*freqstep;
RsignuK1=1/(4*pi)*real((RGfR-RGfA)*(RbathA2-RbathR2)+RGfK*RbathK2)*freqstep;
LsignuR=LsignuR1+conj(LsignuR1(:,end:-1:1));
RsignuR=RsignuR1+conj(RsignuR1(:,end:-1:1));
LsignuK=LsignuK1-conj(LsignuK1(:,end:-1:1));
RsignuK=RsignuK1-conj(RsignuK1(:,end:-1:1));

LsignuA=-conj(LsignuR);
RsignuA=-conj(RsignuR);

for ii=1:length(LInd)
    signuR(LInd(ii),LInd(ii),:)=LsignuR(ii,:);
    signuA(LInd(ii),LInd(ii),:)=LsignuA(ii,:);
    signuK(LInd(ii),LInd(ii),:)=LsignuK(ii,:);
end
for ii=1:length(RInd)
    signuR(RInd(ii),RInd(ii),:)=RsignuR(ii,:);
    signuA(RInd(ii),RInd(ii),:)=RsignuA(ii,:);
    signuK(RInd(ii),RInd(ii),:)=RsignuK(ii,:);
end
end

function n=fermiD(E)
global T
n=1./(exp(E/T)+1);
end