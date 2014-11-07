SP=load('rez_adp.mat');
T=1:size(SP.Pw,2);
Pwf=zak_Pw2+dob_Pw2;

% for i=1:11
%     figure(i)
%      for k=1:3
%         j=(i-1)*3+k;
%         subplot(3,1,k)
%         [hAx,hL1,hL2]=plotyy([T',T',T'],[SP.Pw(j,:)',SP.Pw2(j,:)',Pwf(j,:)'],T,SP.dob_lik2(j,:));
%         ylabel(hAx(1),'Заб. дав., атм.','FontSize',8) % left y-axis
%         ylabel(hAx(2),'Дебит жидк, м^3/сут','FontSize',8) % right y-axis
%         hL2.LineStyle = '--';
%          if k==1
%              legend('Без адаптации Pw','Адаптация Pw','Факт Pw','Дебит факт');
%              set(legend,'Location','west','FontSize',8);
%          end
%         title(SP.WN(k));
%         grid on
%      end;
% end;

  Qo(:,:)=Q(:,3,:);
  Ql(:,:)=Q(:,2,:);
  Qz(:,:)=Q(:,1,:);
  
  Qo2(:,:)=Q2(:,3,:);
  Ql2(:,:)=Q2(:,2,:);
  Qz2(:,:)=Q2(:,1,:);
  
for i=1:11
    figure(i)
     for k=1:3
        j=(i-1)*3+k;
        subplot(3,1,k)
        [hAx,hL1,hL2]=plotyy([T',T',T'],[SP.Pw(j,:)',SP.Pw2(j,:)',Pwf(j,:)'],T,SP.dob_lik2(j,:));
        ylabel(hAx(1),'Заб. дав., атм.','FontSize',8) % left y-axis
        ylabel(hAx(2),'Дебит жидк, м^3/сут','FontSize',8) % right y-axis
        hL2.LineStyle = '--';
         if k==1
             legend('Без адаптации Pw','Адаптация Pw','Факт Pw','Дебит факт');
             set(legend,'Location','west','FontSize',8);
         end
        title(SP.WN(k));
        grid on
     end;
end;