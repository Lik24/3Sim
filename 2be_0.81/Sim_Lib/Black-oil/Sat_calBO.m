function [Zw,Zo,Zg]=Sat_calBO(x2,x3,C,a_s,aw)
%global SwrCUB SorCUB
  
% if isempty(SwrCUB)==1
%     Swr=0.45;
% else
%     Swr=SwrCUB;
% end;
% 
% if isempty(SorCUB)==1
%     Sor=0.10;
% else
%     Sor=SorCUB;
% end;
SwrW=aw(5);
SorW=aw(6);
SorG=aw(7);
SgrG=aw(8);
aw=aw(1:4);
Sg = 1 - x2 - x3;
%   a_s=[2,2];
%  % aw=0.52;
%   aw=0.30;
  
  if C==1
      %Z=interp1(AA(:,1),AA(:,4-B),x2);
      SwW=(x2 - SwrW)./(1 - SorW - SwrW);
      SgG=(Sg - SgrG)./(1 - SgrG - SorG - SwrW);
    
      krw=aw(:,1).*SwW.^a_s(1,1);
      v1w=x2<=(1-SorW);
      v2w=x2>=SwrW;         
      Zw=krw.*v1w.*v2w+aw(:,1).*(v1w==0);
 
      krg=aw(:,4).*SgG.^a_s(1,4);
      v1g=Sg<=(1 - SorG - SwrW);
      v2g=Sg>=SgrG;
      Zg=krg.*v1g.*v2g+aw(:,4).*(v1g==0);
      
      v1ow = 1-x2-SorW>=0;
      v2ow = x2>=SwrW;
      Zow = aw(:,2).*(1-SwW).^a_s(1,2).*v1ow.*v2ow + aw(:,2).*(v2ow==0);
      v1og = 1-Sg-SwrW-SorG>=0;
      v2og = Sg>=SgrG;
      Zog = aw(:,2).*(1-SgG).^a_s(1,3).*v1og.*v2og + aw(:,2).*(v2og==0);
         
      Zo=aw(:,2).*((Zow./aw(:,2) + Zw).*(Zog./aw(:,2) + Zg) - (Zw + Zg));
     
  else
      % Z=interp1(0:0.01:1,D1D(:,B),x2);
      SW=(x2-Swr)./(1-Sor-Swr);
      Y=1./(1-Swr-Sor);
      if B==2
          dkro=-a_s(1,2).*aw(:,2).*Y.*(1-SW).^(a_s(:,2)-1);
          Z=dkro.*(x2<=(1-Sor)).*(x2>Swr);
      else
          dkrw=a_s(1,1)*aw(:,1).*Y.*SW.^(a_s(:,1)-1);
          Z=dkrw.*(x2<(1-Sor)).*(x2>=Swr);
      end;
  end;
  
end


  
