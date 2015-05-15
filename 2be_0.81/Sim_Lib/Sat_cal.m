function [Zw,Zo]=Sat_cal(x2,C,a_s,aw)
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
Swr=aw(5);
Sor=aw(6);
aw=aw(1:2);
%   a_s=[2,2];
%  % aw=0.52;
%   aw=0.30;  
  if C==1
      %Z=interp1(AA(:,1),AA(:,4-B),x2);
      SW=(x2-Swr)./(1-Sor-Swr);
      v1=x2<=(1-Sor);
      v2=x2>=Swr;
     
      krw=aw(:,1).*SW.^a_s(1,1);
      Zw=krw.*v1.*v2+aw(:,1).*(v1==0); 
      kro=aw(:,2).*(1-SW).^a_s(1,2);
      Zo=kro.*v1.*v2+aw(:,2).*(v2==0);    
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


  
