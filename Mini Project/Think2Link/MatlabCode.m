m1 = 10;
m2 = 5;
l1 = 0.2;
l2 = 0.1;
g = 9.81;
zinit = [0, 0.1,0, 0, 0.1, 0];
tspan = [0, 15];
M_fun = @(q2) [(m1 + m2) * l1^2 + m2 * l2 * (l2 + 2 * l1 * cos(q2)), m2 * l2 * (l2 + l1 * cos(q2));
              m2 * l2 * (l2 + l1 * cos(q2)), m2 * (l2^2)];

C_fun = @(q, qdot) [-m2 * l1 * l2 * sin(q(2)) * qdot(2), -m2 * l1 * l2 * sin(q(2)) * (qdot(1) + qdot(2));
                   0, m2 * l1 * l2 * sin(q(2)) * qdot(2)];

G_fun = @(q) [(m1 * l1 * g * cos(q(1))) + (m2 * g * ((l2 * cos(q(1) + q(2))) + (l1 * cos(q(1)))));
             m2 * g * l2 * cos(q(1) + q(2))];
[T, x] = ode45(@(t, z) statespace(t, z, M_fun, C_fun, G_fun), tspan, zinit);


plot(T, x(:, 2), 'r', T, x(:, 5), 'b');
xlabel('Time');
ylabel('Angles q1 and q2');
legend('q1', 'q2');

function y = statespace(t, z, M_fun, C_fun, G_fun)
    y = zeros(size(z));
    qdot=[z(3);z(6)]; 
    q = [z(2);z(5)];  
    ie = [z(1);z(4)];
    q2 = z(5); 
    M = M_fun(q2);
    C = C_fun(q, qdot);
    G = G_fun(q);
    q1_d=0;
    q2_d=0;
    q_d=[q1_d;q2_d];
    M_inverse = inv(M);
       
    Kp = [60;60];%Put Kp,kd and Ki values as you wish.
    Kd = [25;25];
    Ki=[50;50];
    Y = -M_inverse * ((C * qdot) + G) +(Kp .*(q_d-q)) -(Kd .* qdot) +(Ki .* ie);

    
    y(1)=q1_d-z(2);
    y(2)=z(3);
    y(3)=Y(1);
    y(4)=q2_d-z(5);
    y(5)=z(6); 
    y(6)=Y(2);
end




