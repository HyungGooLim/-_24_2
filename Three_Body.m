function dydt = Three_Body(t,u)

global m1 m2 m3

x1= u(1);
y1= u(2);
x1_dot = u(7);
y1_dot = u(8);

x2= u(3);
y2= u(4);
x2_dot = u(9);
y2_dot = u(10);

x3= u(5);
y3= u(6);
x3_dot = u(11);
y3_dot = u(12);

r21 = sqrt((x2-x1)^2 + (y2-y1)^2);
r31 = sqrt((x3-x1)^2 + (y3-y1)^2);
r32 = sqrt((x3-x2)^2 + ((y3-y2)^2));

x1_ddot = m2*(x2-x1)/(r21)^3 + m3*(x3-x1)/(r31)^3;
y1_ddot = m2*(y2-y1)/(r21)^3 + m3*(y3-y1)/(r31)^3;

x2_ddot = m1*(x1-x2)/(r21)^3 + m3*(x3-x2)/(r32)^3;
y2_ddot = m1*(y1-y2)/(r21)^3 + m3*(y3-y2)/(r32)^3;

x3_ddot = m1*(x1-x3)/(r31)^3 + m2*(x2-x3)/(r32)^3;
y3_ddot = m1*(y1-y3)/(r31)^3 + m2*(y2-y3)/(r32)^3;


dydt = [x1_dot; y1_dot;
    x2_dot; y2_dot;
    x3_dot; y3_dot;
    x1_ddot; y1_ddot;
    x2_ddot; y2_ddot;
    x3_ddot; y3_ddot;];

end