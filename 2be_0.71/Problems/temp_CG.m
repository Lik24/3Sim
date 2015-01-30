function [C,G,A2C,A2G,WonC,WonG]=temp_CG(i,C1,G1,A2C1,A2G1,WonC1,WonG1,t)
dy=365*10;
    A2C=A2C1;
    A2G=A2G1;
    
if i<7
    vD=0.1;
elseif i<12
    vD=0.001;
else
    vD=0.00001;
end;

if i==1
    C=C1;
    G=G1;
    WonC=WonC1;
    WonG=WonG1;
elseif i==2
    C=C1*vD;
    G=G1;
    WonC=WonC1*vD;
    WonG=WonG1;
elseif i==3
    C=C1;
    G=G1*vD;
    WonC=WonC1;
    WonG=WonG1*vD;
elseif i==4
    dy=365*10;
    if mod(t,dy)<=0.5*dy
        C=C1*vD;
        G=G1;
        WonC=WonC1*vD;
        WonG=WonG1;
    else
        C=C1;
        G=G1*vD;
        WonC=WonC1;
        WonG=WonG1*vD;
    end;
elseif i==5
    dy=365*1;
    if mod(t,dy)<=0.5*dy
        C=C1*vD;
        G=G1;
        WonC=WonC1*vD;
        WonG=WonG1;
    else
        C=C1;
        G=G1*vD;
        WonC=WonC1;
        WonG=WonG1*vD;
    end;
elseif i==6
    dy=30;
    if mod(t,dy)<=0.5*dy
        C=C1*vD;
        G=G1;
        WonC=WonC1*vD;
        WonG=WonG1;
    else
        C=C1;
        G=G1*vD;
        WonC=WonC1;
        WonG=WonG1*vD;
    end;
elseif i==7
    C=C1*vD;
    G=G1;
    WonC=WonC1*vD;
    WonG=WonG1;
elseif i==8
    C=C1;
    G=G1*vD;
    WonC=WonC1;
    WonG=WonG1*vD;
elseif i==9
    dy=365*10;
    if mod(t,dy)<=0.5*dy
        C=C1*vD;
        G=G1;
        WonC=WonC1*vD;
        WonG=WonG1;
    else
        C=C1;
        G=G1*vD;
        WonC=WonC1;
        WonG=WonG1*vD;
    end;
elseif i==10
    dy=365*1;
    if mod(t,dy)<=0.5*dy
        C=C1*vD;
        G=G1;
        WonC=WonC1*vD;
        WonG=WonG1;
    else
        C=C1;
        G=G1*vD;
        WonC=WonC1;
        WonG=WonG1*vD;
    end;
elseif i==11
    dy=30;
    if mod(t,dy)<=0.5*dy
        C=C1*vD;
        G=G1;
        WonC=WonC1*vD;
        WonG=WonG1;
    else
        C=C1;
        G=G1*vD;
        WonC=WonC1;
        WonG=WonG1*vD;
    end;  
elseif i==12
    C=C1*vD;
    G=G1;
    WonC=WonC1*vD;
    WonG=WonG1;
elseif i==13
    C=C1;
    G=G1*vD;
    WonC=WonC1;
    WonG=WonG1*vD;
elseif i==14
    dy=365*10;
    if mod(t,dy)<=0.5*dy
        C=C1*vD;
        G=G1;
        WonC=WonC1*vD;
        WonG=WonG1;
    else
        C=C1;
        G=G1*vD;
        WonC=WonC1;
        WonG=WonG1*vD;
    end;
elseif i==15
    dy=365*1;
    if mod(t,dy)<=0.5*dy
        C=C1*vD;
        G=G1;
        WonC=WonC1*vD;
        WonG=WonG1;
    else
        C=C1;
        G=G1*vD;
        WonC=WonC1;
        WonG=WonG1*vD;
    end;
else
    dy=30;
    if mod(t,dy)<=0.5*dy
        C=C1*vD;
        G=G1;
        WonC=WonC1*vD;
        WonG=WonG1;
    else
        C=C1;
        G=G1*vD;
        WonC=WonC1;
        WonG=WonG1*vD;
    end; 
end;