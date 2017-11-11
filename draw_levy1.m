format compact

x1 = -5:0.05:5;
x2 = -5:0.05:5;

[X1,X2] = meshgrid(x1,x2);

F = pi/2 * ...
(10*sin(pi*(1+(X1-1)/4)).^2 ...
+ (X1-1).^2/16 .* (1+10*sin(pi*(1+(X2-1)/4)).^2) ...
+ (X2-1).^2/16);

minF = min(min(F));
maxF = max(max(F));

%figure(1), clf
%mesh(X1,X2,F)
%axis tight
%xlabel('x_1')
%ylabel('x_2')
%set(gca, 'XTick',-5:5, 'YTick',-5:5)

figure(1), clf
[C,h] = contour(X1,X2,F, [0:1:9,10:5:50], 'b-');
axis image
xlabel('x_1')
ylabel('x_2')
set(gca, 'XTick',-5:5, 'YTick',-5:5)
clabel(C,h)
hold on
plot(1,1,'r*')
