function [] = save_in_parfor(data,i)
filename = ['Mats',filesep,'pscreen_',num2str(i),'.mat'];
save(filename,'data')

end
% test = rand(10000);
% save('test.mat','test')