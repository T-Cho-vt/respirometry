Params = Ex1_Param_Setup_sim(n,params);

h0 = ones(params.supp,1)/params.supp; % initial guess

% for n=512, true delay is 9.
% Test various delay
Truedelay = params.delay;
DelayTest = [ Truedelay-7 Truedelay Truedelay+7 ];
DelayTest_x = zeros(n,length(DelayTest));
DelayTest_h = zeros(params.supp,length(DelayTest));
for i = 1:length(DelayTest)
    % Different delay
    Params.delay = DelayTest(i);    
    % Run main code
    [blind_x,blind_h,~,blind_res,blind_flag] = Ex1_blind_sim(h0,Params);
%     close all
    DelayTest_x(:,i) = blind_x(:);
    DelayTest_h(:,i) = blind_h(:);
end

Params.hnrm = 1e-6; % Relative error criteria for h
Params.rnrm = 1e-6; % Residual error
%%
figure, hold on
for i = 1:length(DelayTest)
    plot([zeros(DelayTest(i),1);DelayTest_h(:,i)],'linewidth',1.2)
end
hold off

axis([1 120 0 0.08])
legend(['delay=',num2str(DelayTest(1))],['delay=',num2str(DelayTest(2)),' (true)'],['delay=',num2str(DelayTest(3))])
set(gcf,'Position',[300 300 300 200])
