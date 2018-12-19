function [ rVect,vVect ] = coe2rvECI( COE,amu )
%Funzione per passare dai parametri orbitali classici ai vettori velocita'
%e posizione in ECI
a=COE(1);       %Semiasse maggiore [km]
ecc=COE(2);     %Eccentricita'
ainc=COE(3);    %Inclinazione [rad]
gom=COE(4);     %RAAN [rad]
pom=COE(5);     %Argomento del pericentro [rad]
anu=COE(6);     %Anomalia vera [rad]

%Matrici di rotazione
A1=[cos(gom) -sin(gom) 0;
    sin(gom) cos(gom) 0;
    0 0 1];
A2=[1 0 0;
    0 cos(ainc) -sin(ainc);
    0 sin(ainc) cos(ainc)];
A3=[cos(pom) -sin(pom) 0;
    sin(pom) cos(pom) 0;
    0 0 1];
A=A1*A2*A3;

%Terna ECI
g1=[1 0 0];                 %Versore X ECI
g2=[0 1 0];                 %Versore Y ECI
% g3=[0 0 1];                 %Versore Z ECI

%Terna perifocale
p1=A*g1(:);
p2=A*g2(:);
% p3=A*g3(:);

p=a*(1-ecc^2);              %Semilato retto [km]
h=sqrt(p*amu);              %Momento angolare [km^2/s]
r=p/(1+ecc*cos(anu));       %Distanza [km]

%Posizione nel perifocale
rp1=r*cos(anu);
rp2=r*sin(anu);

%Velocita' nel perifocale
vp1=-amu/h*sin(anu);
vp2=amu/h*(ecc+cos(anu));

rVect=rp1*p1+rp2*p2;
vVect=vp1*p1+vp2*p2;
end

