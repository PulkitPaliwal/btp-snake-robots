function out = derivative_novel(t,A,w,i,phi)
diminish = 0.001;
t_d = t(1):diminish:t(length(t));
z = novel_wave(t_d,A,w,i,phi);
out_d = zeros([length(t_d)]);
for i = 1:1:(length(t_d)-1)
    out_d(i,:) = (z(i+1) - z(i))/diminish;
end
out = zeros(length(t));
st = (t(length(t))-t(1))/(length(t)-1);
for i = 1:1:length(t)-1
    out(i) = out_d((1 + (i-1)*st/diminish));
end
out = out(:,1);