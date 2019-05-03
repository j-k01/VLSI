%Transfer Function - VLSI

format long

f = linspace(.01,100,100000);

resultVdb = linspace(1,10000,100000);
resultVp = linspace(1,10000,100000);
s = f.*1j;
%Hz mode
%s= f.2.*pi.*1j;

for z = 1:100000
    
    G1 =1;
    G2 =1;
    L1 = 1.0019*s(z);
    L2 = 1.0019*s(z);
    L3 = 0.2897*s(z);
    C1 = 0.8925*s(z);
    
    T = [G1 0 0 0 1 0 0;
        0 0 0 0 -1 1 1;
        0 0 G2 0 0 -1 0;
        0 0 0 C1 0 0 -1;
        1 -1 0 0 -L1 0 0;
        0 1 -1 0 0 -L2 0;
        0 1 0 -1 0 0 -L3;
        ];
    RHS = [1;0;0;0;0;0;0];
    
    x =T\RHS;
    dbVout = 20*log10(abs(x(3)));
    phaseVout = rad2deg(angle(x(3)));
    resultVdb(z) = dbVout;
    resultVp(z) = phaseVout;


end
    
    hold on
    subplot(2,2,1);
    semilogx(f, resultVdb);
    title('Vout (db)');
    ylabel('Voltage (db)');
    xlabel('Frequency (Hz)');
    subplot(2,2,2);
    semilogx(f, resultVp);
    title('Vout (phase)');
    ylabel('Phase (rad)');
    xlabel('Frequency (Hz)');


%%begin tf calculation

s = [];
D=[];
N=[];
n = 4;
w = exp((2*pi *1i)/(n+1));
for z=0:n
    s(z+1) = w^z;
end
for z=0:n
 
    G1 =1;
    G2 =1;
    L1 = 1.0019*s(z+1);
    L2 = 1.0019*s(z+1);
    L3 = 0.2897*s(z+1);
    C1 = 0.8925*s(z+1);
    
    T = [G1 0 0 0 1 0 0;
        0 0 0 0 -1 1 1;
        0 0 G2 0 0 -1 0;
        0 0 0 C1 0 0 -1;
        1 -1 0 0 -L1 0 0;
        0 1 -1 0 0 -L2 0;
        0 1 0 -1 0 0 -L3;
        ];
    RHS = [1;0;0;0;0;0;0];
    
    x =T\RHS;
   
    D(z+1) = det(T);
    N(z+1) = D(z+1)*x(3);
end



a = [];
for c=0:n
    inter=[];
    for z=0:n 
        inter(z+1) = N(z+1)*(w^(-1*c*z));
    end
    a(c+1) = sum(inter)/(n+1);
end

b = [];
for c=0:n
    inter=[];
    for z=0:n 
        inter(z+1) = D(z+1)*(w^(-c*z));
    end
    b(c+1) = sum(inter)/(n+1);
end


sys = tf(real(flip(a)),real(flip(b)))
%example tf in output and remove zeros
sys = tf([-0.2586 0 -1], [-1.414, -2.306, -2.898, -2])
subplot(2,2,[3,4])
bode(sys)
hold off

