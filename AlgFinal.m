function out = AlgFinal(in)
global a
global b
global c
global seeding 
global stepsize
a = in(1)
b = in(2)
c = in(3)
syms b1 c1 a1
a1 = a
b1 = b
c1 = c
stepsize = in(4)
if stepsize == 0
    stepsize = rand
end
seeding = in(5)
i = 0;
while i < 500
    syms y
    sol1 = vpasolve(-b1 * (y - sin(y)) - c1 * sin(y) == 0, y)
    sol2 = vpasolve(-b1 * (y - sin(y)) - c1 * sin(y) == 0, y, 1)
    sol3 = vpasolve(-b1 * (y - sin(y)) - c1 * sin(y) == 0, y, -1)
    S = [sol1 sol2 sol3]

    BifCheck(S(1,1))
    BifCheck(S(1,2))
    BifCheck(S(1,3))

    i = i + 1
    a1 = a1 + rand
    b1 = b1 + rand
    c1 = c1 + rand
    a = a1 + rand
    b= b1
    c = c1

end
end

function out2 = BifCheck(y)
x = double(y)
iters = 0;
global a
global b
global c
global seeding
global stepsize
while(iters < 2500)
    
    out2 = "";
    iters = iters + 1;
    mat = [((-1 * a) *(cos(x) - 1)) a 0; 1 -1 1; 0 (-1 * b) c];
    
    e = eig(mat);
    v1 = e(1,1);
    v2 = e(2,1);
    v3 = e(3,1);
    
    opposite_signs = false;
    if (real(v3) < 0 & real(v1) > 0) | (real(v3) > 0 & real(v1) < 0)
        opposite_signs = true;
    end
    
    comp_vectors = false;
    if imag(v1) ~= 0
        comp_vectors = true;
    end
     
    if abs(real(v3)) > 0 & imag(v3) == 0 & comp_vectors
        if abs(real(v1)) > 0 & abs(real(v2)) > 0 & real(v1)~= 0 & opposite_signs
              if abs(real(v2)) < abs(v3)
                [V,U] = eig(mat);

                testLinearity = equationstest(V);
                if testLinearity == 0
                    disp("Shilnikov Bifurcation")
                    disp([a, b, c])
                    disp(e)
                    disp(V)
                    out2 = ([a, b, c]);
                    %file output
                    
                    %columns are a b c v1 v2 v3 x0
                    P = [a b c v1 v2 v3 x];
                    G = ["a", "b", "c", "v1", "v2", "v3"];
                    T = table(P);
                    writetable(T, 'shilnikov.csv','WriteMode','append');
                    %check for hopf bifurcations
                    
                    if abs(real(v1)) < 0.05 & abs(real(v2)) < 0.05
                        out2 = hopf(x);
                        
                    end
                end
             end
        end
    end
    
    
    
    %linear search
    %a = a + 0.05;
    %b = b + 0.05;
    %c = c + 0.05;
    
    if seeding == 0 
        a = a + stepsize;
        b = b + stepsize;
        c = c + stepsize;
    else
        a = a + (stepsize * rand);
        b = b + (stepsize * rand);
        c = c + (stepsize * rand);
    
    end
    end
end 

function out3 = equationstest(V)
u1 = V(1,1);
u2 = V(2, 1);
u3 = V(3,1);

v1 = V(1,2);
v2 = V(2, 2);
v3 = V(3,2);

w1 = V(1,3);
w2 = V(2, 3);
w3 = V(3,3);


syms f g h 
sol = solve([f*u1 - g*v1 + h * w1 == 0, f*u2 + g*v2 - w2*h == 0, u3 + v3 + u3 == 0], [f,g,h]);
solA = sol.f;
solB = sol.g;
solC = sol.h;
if solA == 0 & solA == solB & solB == solC
    out3 = 0;
else
    out3 = 1;

end
end

function hop = hopf(x)
    global a
    global b
    global c
    hop = 0;
    a1 = a;
    b1 = b;
    c1 = c;
    mat = [0 a1 0; 1 -1 1; 0 -b1 c1]

    e = eig(mat);
    v1 = e(1,1);
    v2 = e(2,1);
    v3 = e(3,1);
    pastv = 0 + 0i;
    pastv2 = 0 + 0i;
    pastv3 = 0;
    while real(v1) ~= 0
        %a = a - (rand * rand);
        b1 = b1 + (rand * rand);
        c1 = c1 + (rand * rand);
        mat = [0 a1 0; 1 -1 1; 0 -b1 c1];
        e = eig(mat);
        v1 = e(1,1);
        v2 = e(2,1);
        v3 = e(3,1);

        if real(v1) == 0
            if real(v2) == 0 
                disp("Hopf Bifurcation")
                disp([a1 b1 c1])
                disp ([v1 v2 v3])
                hop = ([a1 b1 c1])
                
                break
            end
        end
        
        if real(v1) > 0
            if real(pastv < 0)
                disp("Hopf Bifurcation Near")
                disp([a1 b1 c1])
                disp("Sign changed from - to +")
                disp ([v1 v2 v3])
                disp([pastv pastv2])
                P = [a b c v1 v2 v3 pastv pastv2 pastv3 x real(0)];
                T = table(P);
                writetable(T, 'hopf.csv','WriteMode','append');
                hop = ([a1 b1 c1])
                break
            end
        end
        
        if real(v1) < 0
            if real(pastv) > 0
                disp("Hopf Bifurcation Near")
                disp([a1 b1 c1])
                disp("Sign changed from + to -")
                disp ([v1 v2 v3])
                disp ([pastv pastv2])
                %1 is + -> -
                P = [a b c v1 v2 v3 pastv pastv2 pastv3 x real(1)];
                T = table(P);
                writetable(T, 'hopf.csv','WriteMode','append');
                hop = ([a1 b1 c1])
                break
            end
        end
        
        if (real(v3) < 0 & real(v1) > 0) | (real(v3) > 0 & real(v1) < 0)
            opposite_signs = true;
        end
        
        comp_vectors = false;
    if imag(v1) ~= 0
        comp_vectors = true;
    end
     
    if abs(real(v3)) > 0 & imag(v3) == 0 & comp_vectors
        if abs(real(v1)) > 0 & abs(real(v2)) > 0 & real(v1)~= 0 & opposite_signs
              if abs(real(v2)) < abs(v3)
                [V,U] = eig(mat);

                testLinearity = equationstest(V);
                if testLinearity == 0
                    disp("Shilnikov Bifurcation")
                    disp([a, b, c])
                    disp(e)
                    disp(V)
                    out2 = ([a, b, c]);
                    %file output
                    
                    %columns are a b c v1 v2 v3 x0
                    P = [a b c v1 v2 v3 x];
                    G = ["a", "b", "c", "v1", "v2", "v3"];
                    T = table(P);
                    writetable(T, 'shilnikov.csv','WriteMode','append');
                    %check for hopf bifurcations
                end
              end
        end
    end

        pastv = v1;
        pastv2 = v2;
        pastv3 = v3;
        
    
    end
end 