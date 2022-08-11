clear;clc;close all;
u = 0:1/100:2;
t = 1:length(u);
t = t';
y_m = zeros(length(u),1);

for i = 1:length(u)
    u_l = ones(size(t,1),1)*u(i);
    y = odezva_2021(192291,u_l,t);
    y_m(i) = mean(y(end-5:end));
end

plot(u,y_m)

%% vyjde do cca 0.9