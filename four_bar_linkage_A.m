a= 1;
b= 5.5;
c= 13.642;
d= 13.792;

K1= d/a;
K2= d/c;
K3= (a^2 - b^2 + c^2 + d^2)/(2*a*c);
K4= d/b;
K5= (c^2 - d^2 - a^2 - b^2)/(2*a*b);

theta2 = linspace(0, 2*pi,360);
theta3 = zeros(size(theta2));
theta4 = zeros(size(theta2));
trans_angle = zeros(size(theta2));

%% Finding Theta3 and Theta4 for any angle Theta2
for pos = 1:360
    i = theta2(pos);
    A = cos(i) - K1 - K2*cos(i) + K3;
    B = -2*sin(i);
    C= K1 - (K2+1)*cos(i) + K3;
    D= cos(i) - K1 + K4*cos(i) + K5;
    E= -2*sin(i);
    F= K1 + (K4-1)*cos(i) + K5;

    theta4(pos) = 2*atan((-B-sqrt(B^2-4*A*C))/(2*A));
    theta3(pos) = 2*atan((-E-sqrt(E^2-4*D*F))/(2*D));
    trans_angle(pos) = abs(theta4(pos) - theta3(pos));
end

%% Given Input & Output force vectors and induced Ang_vel2 calculate other Ang_vel
w2=.1;
rin=14.25;
rout=1.8;

w3=zeros(size(theta2));
w4=zeros(size(theta2));
mA=zeros(size(theta2));     % remember that mechanical advantage trends to infinity when our pressure angle changes sign.

for pos = 1:360
    w3(pos)=(a*w2*sin(theta4(pos)-theta2(pos)))/(b*sin(theta3(pos)-theta4(pos)));
    w4(pos)=(a*w2*sin(theta2(pos)-theta3(pos)))/(c*sin(theta4(pos)-theta3(pos)));
    mA(pos)=abs((c*sin(theta4(pos)-theta3(pos))*rin)/ ...
        (a*sin(theta2(pos)-theta3(pos))*rout) ...
        );
end

%% Given Center of Gravity lengths & andgles and input accaleration, calculate angular acceleration
% and acceleration of each links' center of gravity
alph2=24.12;
s=6.625;
p=2.75;
u=6.204;
delt2=180*pi/180;
delt3=0;
delt4=0;

alph3=zeros(size(theta2));
alph4=zeros(size(theta2));
A_G2_x=zeros(size(theta2));
A_G2_y=zeros(size(theta2));
A_G3_x=zeros(size(theta2));
A_G3_y=zeros(size(theta2));
A_G4_x=zeros(size(theta2));
A_G4_y=zeros(size(theta2));

for pos = 1:360
    A=c*sin(theta4(pos));
    B=b*sin(theta3(pos));
    C=a*alph2*sin(theta2(pos))+a*((w2)^2)*cos(theta2(pos))+b*((w3(pos))^2)*cos(theta3(pos))-...
        c*(w4(pos))^2*cos(theta4(pos));
    D=c*cos(theta4(pos));
    E=b*cos(theta3(pos));
    F=a*alph2*cos(theta2(pos))-a*w2^2*sin(theta2(pos))-b*(w3(pos))^2*sin(theta3(pos))+...
        c*(w4(pos))^2*sin(theta4(pos));
    alph3(pos)=(C*D-A*F)/(A*E-B*D);
    alph4(pos)=(C*E-B*F)/(A*E-B*D);

    A_G2_x(pos)=-(s*alph2*sin(theta2(pos)+delt2)+s*w2^2*cos(theta2(pos)+delt2));
    A_G2_y(pos)=s*alph2*cos(theta2(pos)+delt2)-s*w2^2*sin(theta2(pos)+delt2);
    
    A_G3_x(pos)=-(a*alph2*sin(theta2(pos))+a*w2^2*cos(theta2(pos)))...
        -(p*alph3(pos)*sin(theta3(pos)+delt3)+p*w3(pos)^2*cos(theta3(pos)+delt3));
    A_G3_y(pos)=(a*alph2*cos(theta2(pos))-a*w2^2*sin(theta2(pos))+...
        p*alph3(pos)*cos(theta3(pos)+delt3)-p*w3(pos)^2*sin(theta3(pos)+delt3));

    A_G4_x(pos)=-(u*alph4(pos)*sin(theta4(pos)+delt4)+u*w4(pos)^2*cos(theta4(pos)+delt4));
    A_G4_y(pos)=u*alph4(pos)*cos(theta4(pos)+delt4)-u*w4(pos)^2*sin(theta4(pos)+delt4);
end

%% Given mass of each link, determine moment of inertia. Additionally, find position vectors, force vectors and torque vectors
%masses {blobs}
m2=7.513*(10^-4);
m3=2.332*(10^-5);
m4=7.642*(10^-4);

%I-values {blobs*in^2}
I_G2=0.01601;
I_G3=5.88*(10^-5);
I_G4=0.01637;

%R-values
R_12_x=zeros(size(theta2));
R_12_y=zeros(size(theta2));

R_32_x=zeros(size(theta2));
R_32_y=zeros(size(theta2));

R_23_x=zeros(size(theta2));
R_23_y=zeros(size(theta2));

R_43_x=zeros(size(theta2));
R_43_y=zeros(size(theta2));

R_14_x=zeros(size(theta2));
R_14_y=zeros(size(theta2));

R_34_x=zeros(size(theta2));
R_34_y=zeros(size(theta2));

%populating R-values
for pos = 1:360
        %link 2
    R_12_x(pos)=s*cos(theta2(pos)+delt2+pi);
    R_12_y(pos)=s*sin(theta2(pos)+delt2+pi);
    R_32_x(pos)=a*cos(theta2(pos))+R_12_x(pos);
    R_32_y(pos)=a*sin(theta2(pos))+R_12_y(pos);
        %link 3
    R_23_x(pos)=p*cos(theta3(pos)+delt3+pi);
    R_23_y(pos)=p*sin(theta3(pos)+delt3+pi);
    R_43_x(pos)=b*cos(theta3(pos))+R_23_x(pos);
    R_43_y(pos)=b*sin(theta3(pos))+R_23_y(pos);
        %link 4 
    R_14_x(pos)=u*cos(theta4(pos)+delt4+pi);
    R_14_y(pos)=u*sin(theta4(pos)+delt4+pi);
    R_34_x(pos)=c*cos(theta4(pos))+R_14_x(pos);
    R_34_y(pos)=c*sin(theta4(pos))+R_14_y(pos);
end

F_12_x=zeros(size(theta2));
F_12_y=zeros(size(theta2));
F_32_x=zeros(size(theta2));
F_32_y=zeros(size(theta2));
F_43_x=zeros(size(theta2));
F_43_y=zeros(size(theta2));
F_14_x=zeros(size(theta2));
F_14_y=zeros(size(theta2));
T_12=zeros(size(theta2));
Mag_12=zeros(size(theta2));
Mag_32=zeros(size(theta2));
Mag_34=zeros(size(theta2));
Mag_14=zeros(size(theta2));


for pos = 1:360
    N=[ 1 0 1 0 0 0 0 0 0; 
        0 1 0 1 0 0 0 0 0; 
        -R_12_y(pos) -R_12_x(pos) -R_32_y(pos) -R_32_x(pos) 0 0 0 0 1;
        0 0 -1  0 1 0 0 0 0;
        0 0 0 -1 0 1 0 0 0;
        0 0 R_23_y(pos) -R_32_x(pos) -R_43_y(pos) R_43_x(pos) 0 0 0;
        0 0 0 0 -1 0 1 0 0;
        0 0 0 0 0 -1 0 1 0;
        0 0 0 0 R_34_y(pos) -R_34_x(pos) -R_14_y(pos) R_14_x(pos) 0];

    M=[ m2*A_G2_x(pos);
        m2*A_G2_y(pos);
        I_G2*alph2;
        m3*A_G3_x(pos);
        m3*A_G3_y(pos);
        I_G3*alph3(pos);
        m4*A_G4_x(pos)+103*cos(pi-theta4(pos));
        m4*A_G4_y(pos)+103*sin(pi-theta4(pos));
        I_G4*alph4(pos)-16*cos(pi-theta4(pos))*103*sin(pi-theta4(pos))];

    Force_vec=N\M;
    Force_vec=Force_vec';
    F_12_x(pos)=Force_vec(1);
    F_12_y(pos)=Force_vec(2);
    Mag_12(pos)=sqrt(F_12_x(pos)^2+F_12_y(pos)^2);

    F_32_x(pos)=Force_vec(3);
    F_32_y(pos)=Force_vec(4);
    Mag_32(pos)=sqrt(F_32_x(pos)^2+F_32_y(pos)^2);

    F_43_x(pos)=Force_vec(5);
    F_43_y(pos)=Force_vec(6);
    Mag_43(pos)=sqrt(F_43_x(pos)^2+F_43_y(pos)^2);

    F_14_x(pos)=Force_vec(7);
    F_14_y(pos)=Force_vec(8);
    Mag_14(pos)=sqrt(F_14_x(pos)^2+F_14_y(pos)^2);

    T_12(pos)=Force_vec(9);
end


%% poltting results. user can pick and choose!
while true
        
    disp('What behavior would you like to plot?');
    disp('1: position');
    disp('2: velocity');
    disp('3: acceleration');
    disp('4: mechanical adantage');
    disp('5: pin force reactions')
    choice= input('Enter 1-5. (press q to quit)',"s");
    disp('----------------------------------------------------');

    if choice == 'q'
        disp('exiting...')
        break
    end

    switch choice
        case "1"      %position, theta_3 and theta_4
            angle_choice=input('Do you want theta 3, theta 4, or transmission angle? (3, 4, t)', "s");
            
            switch angle_choice
                case "3"
                    plot(theta2, theta3);
                    xlabel('Theta2 (rad)');
                    ylabel('Theta3 (rad)');
                    title('Angular Position of Link 3')
    
                case "4"
                    plot(theta2, theta4);
                    xlabel('Theta2 (rad)');
                    ylabel('Theta4 (rad)');
                    title('Angular Position of Link 4')

                case "t" 
                    plot(theta2, trans_angle);
                    xlabel('Theta2 (rad)');
                    ylabel('Transmission angle (rad)');
                    title('Direction of Output Force');

                otherwise
                    disp('');
                    disp('!_____!');
                    disp('unknown input. going back to top.');
                    disp(' ');
                    continue
            end
    
        case "2"    %velocities
            vel_choice=input('Do you want w3 or w4? (3 or 4)', "s");
            
            switch vel_choice
                case "3"
                    plot(theta2, w3);
                    xlabel('Theta2 (rad)');
                    ylabel('Omega3 (rad/s)');
                    title('Angular Velocity of Link 3');
                case "4"
                    plot(theta2, w4);
                    xlabel('Theta2 (rad)');
                    ylabel('Omega4 (rad/s)');
                    title('Angular Velocity of Link 4');
                otherwise 
                    disp('');
                    disp('!_____!');
                    disp('unknown input. going back to top.');
                    disp(' ');
                    continue
            end
        
        case "3"      %acceleration, angular and center of mass
            disp('Which link do you want to know about?');
            accel_choice=input('Link 2, Link 3, or Link 4 (2, 3, 4)', "s");

            switch accel_choice
                case "2"
                    plot(theta2, A_G2_x, theta2, A_G2_y);
                    xlabel('Theta2 (rad)');
                    ylabel('Linear Acceleration (in/s^2)');
                    legend('x-component', 'y-component');
                    title('Linear Acceleration of Link 2');

                case "3" 
                    ang_or_lin=input('Angular or Linear? (a/l)', 's');
                    if ang_or_lin == "a"
                        plot(theta2, alph3);
                        xlabel('Theta2 (rad)');
                        ylabel('Alpha (rad/s^2)');
                        title('Angular Acceleration of Link 3');
                    elseif ang_or_lin == "l"
                        plot(theta2, A_G3_x, theta2, A_G3_y);
                        xlabel('Theta2 (rad)');
                        ylabel('Linear Acceleration (in/s^2)');
                        legend('x-component', 'y-component');
                        title('Linear Acceleration of Link 3');
                    end

                case "4"
                    ang_or_lin=input('Angular or Linear? (a/l)', 's');
                    if ang_or_lin == "a"
                        plot(theta2, alph4);
                        xlabel('Theta2 (rad)');
                        ylabel('Alpha (rad/s^2)')
                        title('Angular Acceleration of Link 4');
                    elseif ang_or_lin == "l"
                        plot(theta2, A_G4_x, theta2, A_G4_y);
                        xlabel('Theta2 (rad)');
                        ylabel('Linear Acceleration (in/s^2)');
                        legend('x-component', 'y-component');
                        title('Linear Acceleration of Link 4');
                    end
                    
                otherwise
                    disp('');
                    disp('!_____!');
                    disp('unknown input. going back to top.');
                    disp(' ');
                    continue
            end

        case "4"
            plot(theta2, mA);
            xlabel('Theta2 (rad)');
            ylabel('Mechanical Advantage');
            ylim([-10 250]);        % limit the range of the window so the peaks don't skew our scale
    
        case "5"
            choice_pin=input('which pin piques your interest? (1-4)', "s");
            switch choice_pin
                case "1"
                    plot(theta2, F_12_x, theta2, F_12_y, theta2, Mag_12);
                    xlabel('Theta2 (rad)');
                    ylabel('Forces on pin (lbf)');
                    legend('x-component', 'y-component', 'Magnitude');
                    title('Reaction forces on pin 1');
                    
                case "2"
                    plot(theta2, F_32_x, theta2, F_32_y, theta2, Mag_32);
                    xlabel('Theta2 (rad)');
                    ylabel('Forces on pin (lbf)');
                    legend('x-component', 'y-component', 'Magnitude');
                    title('Reaction forces on pin 2');
                    
                case "3"
                    plot(theta2, F_43_x, theta2, F_43_y, theta2, Mag_43);
                    xlabel('Theta2 (rad)');
                    ylabel('Forces on pin (lbf)');
                    legend('x-component', 'y-component', 'Magnitude');
                    title('Reaction forces on pin 3');
                    
                case "4"
                    plot(theta2, F_14_x, theta2, F_14_y, theta2, Mag_14);
                    xlabel('Theta2 (rad)');
                    ylabel('Forces on pin (lbf)');
                    legend('x-component', 'y-component', 'Magnitude');
                    title('Reaction forces on pin 4');
                    
                otherwise
                    disp('');
                    disp('!_____!');
                    disp('unknown input. going back to top.');
                    disp(' ');
                    continue
            end
        otherwise
            disp('');
            disp('!_____!');
            disp('unknown input. going back to top.');
            disp(' ');
            continue
    end
end