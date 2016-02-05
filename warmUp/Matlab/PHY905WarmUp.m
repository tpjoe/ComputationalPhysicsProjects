h = 0.0000001:0.0000001:1;
result = zeros(length(h),1);
for i = 1:length(h)
    result(i) = daTan(sqrt(2),h(i),1,'single');
end
error = log10(3*abs(result-1/3));
[Y,I] = min(error);
hold on
plot(log10(h),error);
%plotyy(h,result,h,error)
%legend('Results','Error');

result = zeros(length(h),1);
for i = 1:length(h)
    result(i) = daTan(sqrt(2),h(i),1,'double');
end
error = log10(3*abs(result-1/3));
[Y,I] = min(error);
hold on
plot(log10(h),error);

result = zeros(length(h),1);
for i = 1:length(h)
    result(i) = daTan(sqrt(2),h(i),2,'single');
end
error = log10(3*abs(result-1/3));
[Y,I] = min(error);
hold on
plot(log10(h),error);

result = zeros(length(h),1);
for i = 1:length(h)
    result(i) = daTan(sqrt(2),h(i),2,'double');
end
error = log10(3*abs(result-1/3));
[Y,I] = min(error);
hold on
plot(log10(h),error);

legend('1,single','1,double','2,single','2,double');