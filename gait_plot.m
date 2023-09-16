t = 0:0.01:2;
w = 2*pi;
A = 2;
generate_plot_for_delta_variation = true;
generate_plot_for_link_angle_variation = false;
if generate_plot_for_link_angle_variation == true
    delta = deg2rad(3);
    Z = zeros([9 length(t)]);
    for i = 1:1:9
        Z(i,:) = novel_wave(t,A,w,i,delta);
    end
    for i = 1:1:9
        plot(t,Z(i,:));
        grid on;
        hold on;
    end
    legend('i=1', 'i=2','i=3','i=4','i=5','i=6','i=7','i=8','i=9');
    title("Variation of Link Angles")
    hold on;
elseif generate_plot_for_delta_variation == true
   delta = [2,8,16,32];
   delta = deg2rad(delta) ;
   Z = zeros([length(delta) length(t)]);
   for i = 1:1:length(delta)
         Z(i,:) = novel_wave(t,A,w,2,delta(i));
         plot(t,Z(i,:));
         grid on;
         hold on;
   end
   legend('\delta = 2','\delta = 8','\delta = 16','\delta = 32');
   title("Variation of Second Link angle with \delta")
end