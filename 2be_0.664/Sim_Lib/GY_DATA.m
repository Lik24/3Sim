function [GY_Data]=GY_DATA(bnd_XY,bnd_Z)
nz=size(bnd_Z);
nxy=size(bnd_XY);

GY_Kz=1.0*ones(nz)*8.64;
GY_Pz=100*ones(nz);
GY_Swz=1*ones(nz);

GY_Kz(bnd_Z~=2)=0;
GY_Kz(bnd_Z~=1)=0;

GY_Kxy=0*ones(nxy)*8.64;
GY_Pxy=100*ones(nxy);
GY_Swxy=1*ones(nxy);

GY_Data.GY_Kz=GY_Kz;
GY_Data.GY_Pz=GY_Pz;
GY_Data.GY_Swz=GY_Swz;

GY_Data.GY_Kxy=GY_Kxy;
GY_Data.GY_Pxy=GY_Pxy;
GY_Data.GY_Swxy=GY_Swxy;