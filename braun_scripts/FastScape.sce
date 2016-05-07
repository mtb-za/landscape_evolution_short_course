global donor ndon stack nstack;

function find_stack(ij)
global donor ndon stack nstack;
    for k=1:ndon(ij)
        ijk=donor(ij,k);
        nstack=nstack+1;
        stack(nstack)=ijk;
        find_stack(ijk)
    end
endfunction

nx=51;
ny=51;
nn=nx*ny;
h=rand(nn,1);
hplot=zeros(nx,ny);
rec=zeros(nn,1);
len=zeros(nn,1);
ndon=zeros(nn,1);
donor=zeros(nn,8);
stack=zeros(nn,1);
nstack=0;
a=zeros(nn,1);

xl=100000;
yl=100000;
dx=xl/(nx-1);
dy=yl/(ny-1);
dt=10000;
nstep=100;
kf=0.0001;
u=0.002;
n=1;
m=0.4;
precipitation=1;
kd=0.1;

for istep=1:nstep

    for ij=1:nn
        rec(ij)=ij;
        len(ij)=0;
    end

    for j=2:ny-1
        for i=2:nx-1
            ij=i+(j-1)*nx;
            smax=0.000001;
            for jj=-1:1
                for ii=-1:1
                    iii=i+ii;
                    jjj=j+jj;
                    ijk=iii+(jjj-1)*nx;
                    if ijk~=ij
                        ll=sqrt((dx*ii)^2+(dy*jj)^2);
                        slope=(h(ij)-h(ijk))/ll;
                        if slope>smax
                            smax=slope;
                            rec(ij)=ijk;
                            len(ij)=ll;
                        end
                    end
                end
            end
        end
    end

    for ij=1:nn
        ndon(ij)=0;
    end

    for ij=1:nn
        if rec(ij)~=ij
            ijk=rec(ij);
            ndon(ijk)=ndon(ijk)+1;
            donor(ijk,ndon(ijk))=ij;
        end
    end

    nstack=0;
    for ij=1:nn
        if rec(ij)==ij
            nstack=nstack+1;
            stack(nstack)=ij;
            find_stack(ij);
        end
    end

    for ij=1:nn
        a(ij)=precipitation*dx*dy;
    end

    for ij=nn:-1:1
        ijk=stack(ij);
        if rec(ijk)~=ijk
            a(rec(ijk))=a(rec(ijk))+a(ijk);
        end
    end

    for j=2:ny-1
        for i=2:nx-1
            ij=i+(j-1)*nx;
            h(ij)=h(ij)+u*dt;
        end
    end

    for ij=1:nn
        ijk=stack(ij);
        ijr=rec(ijk);
        if ijr~=ijk
            ll=len(ijk);
            fact=kf*dt*a(ijk)^m/ll;
            h(ijk)=(h(ijk)+fact*h(ijr))/(1+fact);
        end
    end

    factx=kd*dt/dx**2;
    facty=kd*dt/dy**2;
    for j=2:ny-1
        for i=2:nx-1
           ij=i+(j-1)*nx;
           imj=i-1+(j-1)*nx;
           ipj=i+1+(j-1)*nx;
           ijm=i+(j-2)*nx;
           ijp=i+j*nx;
           h(ij)=h(ij)+factx*(h(ipj)-2*h(ij)+h(imj))+facty*(h(ijp)-2*h(ij)+h(ijm))
        end
    end

    for j=1:ny
        for i=1:nx
            ij=i+(j-1)*nx;
            hplot(i,j)=h(ij);
        end
    end

    x=linspace(0,nx-1,nx);
    y=linspace(0,ny-1,ny);
    clf;
    contourf(x,y,hplot,16);
    contour(x,y,hplot,16);
    disp(istep);

end
