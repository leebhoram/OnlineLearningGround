function q = quatprod(q1, q2)

q = [q1(1) -q1(2) -q1(3) -q1(4); 
     q1(2)  q1(1) -q1(4)  q1(3); 
     q1(3)  q1(4)  q1(1) -q1(2); 
     q1(4) -q1(3)  q1(2)  q1(1)]*q2;

end

