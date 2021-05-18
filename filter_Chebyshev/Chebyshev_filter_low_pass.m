%CHEBYSHEV'S APPROXIMATION FILTER TYPE I - LOW PASS
%The params are:
%f_c - Cutoff Frecuency
%ripple_db - Ripple Db'S
%attenuation_db - Attenuation In The Eliminated Band 
%f_r -Frecuency Attenuation In The Eliminated Band
%Internal it has 2 features more, the resistor size and language for display
%Author Jorge Benavides, rhyloot@gmail.com
%2021-05-17
function filter_low_pass(f_c,ripple_db,attenuation_db,f_r);
%% Introduce The Basics Params
% f_c = 1000;             % Cutoff Frecuency
% ripple_db = 2;      % Ripple Db'S
% attenuation_db = 40; % Attenuation In The Eliminated Band 
% f_r = 1300;             % Frecuency Attenuation In The Eliminated Band
% R = 1e4;                %Design Resistor R
% R1 = 1e4;               %Design Resistor R
%%Set Lang
%Options Available = [Spanish English]
%For Publishing Purposes I Recommend Use English Some Tools Like LaTeX Fails With Accents.
lang = "english";

%%Ripple Calc
syms epsilon
ripple = solve(10^(ripple_db/20) == sqrt(1+epsilon^2));
ripple = ripple(find(ripple>0));

if lang == "spanish"
    fprintf('El rizado es: %f\n',double(ripple))
elseif lang == "english"
    fprintf('The ripple is: %f\n',double(ripple))
end

%%Attenuation Calc
attenuation = 10^(attenuation_db/20);
if lang == "spanish"
    fprintf('Atenuación: %f\n',double(attenuation))
elseif lang == "english"
    fprintf('Attenuation: %f\n',double(attenuation))
end

%%Filter Order
g = double(sqrt((attenuation^2-1)/(ripple^2)));
num = log(g+sqrt(g^2-1));
dem = log(((2*pi*f_r)/(2*pi*f_c))+sqrt((((2*pi*f_r)/(2*pi*f_c))^2-1)));
real_filter_order = num/dem;
filter_order = ceil(num/dem);

if lang == "spanish"
    fprintf('Orden del filtro: %f\n',real_filter_order)
    fprintf('Aproximación del filtro: %f\n',filter_order)
elseif lang == "english"
    fprintf('Real filter order: %f\n',real_filter_order)
    fprintf('Filter order approximation: %f\n',filter_order)
end


%%Coefficients Calc
if lang == "spanish"
    fprintf('Coeficientes:\n')
elseif lang == "english"
    fprintf('Coefficients:\n')
end
%u coeficient
mu = double(ripple^(-1)+sqrt(ripple^(-2)+1));
if lang == "spanish"
    fprintf('Coeficiente u: %f\n',mu)
elseif lang == "english"
    fprintf('Coeficient u: %f\n',mu)
end

%minor axis
a = (mu^(1/filter_order)-mu^(-1/filter_order))/2;
if lang == "spanish"
    fprintf('Coeficiente a: %f\n',a)
elseif lang == "english"
    fprintf('Coeficient a: %f\n',a)
end

%major axis
b = (mu^(1/filter_order)+mu^(-1/filter_order))/2;
if lang == "spanish"
    fprintf('Coeficiente b: %f\n',b)
elseif lang == "english"
    fprintf('Coeficient b: %f\n',b)
end

%%Poles Calcs
list_poles = [];
if lang == "spanish"
    fprintf('Los polos son:\n')
elseif lang == "english"
    fprintf('The poles are:\n')
end
for k = 0:filter_order-1
    pole = -a*sin((pi/(2*filter_order))*(1+2*k))+i*b*cos((pi/(2*filter_order))*(1+2*k));
    if round(imag(pole),10) == 0
        pole = -1*real(pole);
    end
    list_poles = [list_poles;pole];
    if lang == "spanish"
        if imag(list_poles(k+1)) < 0
            fprintf('Polo %d: %f%fi\n',k+1,real(list_poles(k+1)),imag(list_poles(k+1)))
        else
            fprintf('Polo %d: %f+%fi\n',k+1,real(list_poles(k+1)),imag(list_poles(k+1)))
        end
    elseif lang == "english"
        if imag(list_poles(k+1)) < 0
            fprintf('Pole %d: %f%fi\n',k+1,real(list_poles(k+1)),imag(list_poles(k+1)))
        else
            fprintf('Pole %d: %f+%fi\n',k+1,real(list_poles(k+1)),imag(list_poles(k+1)))
        end
    end
end

%%Pair Of Poles Calcs
if lang == "spanish"
    fprintf('Pares de polos, frecuencias, calidades:\n')
elseif lang == "english"
    fprintf('Pairs of poles, frequencies, qualities:\n')
end
%For even order
     pairs = [];
     frequencies = [];
     qualities = [];
if rem(filter_order,2) == 0
     for k = 1:filter_order/2
         pair = [1 -1*2*real(list_poles(k)) real(list_poles(k))^2+imag(list_poles(k))^2];
         pairs = [pairs;pair];
         frequency = sqrt(real(list_poles(k))^2+imag(list_poles(k))^2);
         frequencies = [frequencies;frequency];
         quality = frequencies(k)/(-1*2*(real(list_poles(k))));
         qualities = [qualities; quality];
         fun = printout_TF(tf([pairs(k,:)],[1]));
         fprintf("%s --> omega %d = %f --> quality %d = %f\n",fun,k,frequencies(k),k,qualities(k))
     end
     
else
%For odd order
for k = 1:(filter_order/2)+0.5
    if isreal(list_poles(k))
        pair = [0 1 list_poles(k)];
    else
        pair = [1 -1*2*real(list_poles(k)) real(list_poles(k))^2+imag(list_poles(k))^2];
    end
         pairs = [pairs;pair];
         frequency = sqrt(real(list_poles(k))^2+imag(list_poles(k))^2);
         frequencies = [frequencies;frequency];
         quality = frequencies(k)/(-1*2*(real(list_poles(k))));
         qualities = [qualities; quality];
         fun = printout_TF(tf([pairs(k,:)],[1]));
         fprintf("%s --> omega %d = %f --> quality %d = %f\n",fun,k,frequencies(k),k,qualities(k))
     end
end

%%Gain Of Transfer Function
if rem(filter_order,2) == 0
    H_j0 = double(sqrt(1/(1+ripple^2)));
    H_j0
    for k = 1:filter_order/2
    H_j0 = H_j0 * (frequencies(k))^2;
    end
    H_0 = H_j0;
else
    H_j0 = 1;
    for k = 1:(filter_order/2)+1
    H_j0 = H_j0 * (frequencies(k))^2;
    end
    H_0 = H_j0;
end

if lang == "spanish"
    fprintf('Ganancia: %f\n',H_0)
elseif lang == "english"
    fprintf('Gain: %f\n',H_0)
end

%%Normalized Transfer Function
sys = 1;
if rem(filter_order,2) == 0
    for k = 1:(filter_order/2)
        sys = sys*tf([1],[1 frequencies(k)/qualities(k) frequencies(k)^2]);
    end
    sys = sys*H_0;
else
    for k = 1:(filter_order/2)+1
        if (k == (filter_order/2)+0.5)
            sys = sys*tf([1],[1 frequencies(k)]);    
        else
            sys = sys*tf([1],[1 frequencies(k)/qualities(k) frequencies(k)^2]);
        end
    end
    sys = sys*H_0;
end

if lang == "spanish"
    fprintf('Función de transferencia normalizada:')
    zpk(sys)
elseif lang == "english"
    fprintf('Normalized transfer function:')
    zpk(sys)
end

%%Transfer Function
w_c = 2*pi*f_c;
sys = 1;
if rem(filter_order,2) == 0
    for k = 1:(filter_order/2)
        sys = sys*tf([1],[1  ((frequencies(k)*w_c)/(qualities(k))) (frequencies(k)^2*w_c^2)]);
    end
    sys = sys*H_0*w_c^2;
else
    for k = 1:(filter_order/2)+1
        if (k == (filter_order/2)+0.5)
            sys = sys*tf([w_c^2],[1 frequencies(k)*w_c]);    
        else
            sys = sys*tf([1],[1  ((frequencies(k)*w_c)/(qualities(k))) (frequencies(k)^2*w_c^2)]);
        end
    end
    sys = sys*H_0;
end
%Revisar este trozo
if lang == "spanish"
    fprintf('Función de transferencia')
    zpk(sys)
elseif lang == "english"
    fprintf('Transfer function:')
    zpk(sys)
end

%%Values of resistors and capacitors for implementation
[Z,P,K] = zpkdata(sys,'v');
d = 0;
if rem(filter_order,2) == 0
    for k = 1:2:(filter_order)
        d = d + 1;
        polynomial = poly([P(k) P(k+1)]);
        C2 = 2/(R*polynomial(2));
        C5 = 1/(R^2*C2*polynomial(3));
        if lang == "spanish"
            fprintf('\n Par %d:\n Condensador 2: %e\n Condensador 5: %e\n',d,C2,C5)
        elseif lang == "english"
            fprintf('\n Par %d:\n Capacitor 2: %e\n Capacitor 5: %e\n',d,C2,C5)
        end 
    end
else
    ugly_pole = P(find(imag(P)==0));
    P(find(imag(P) == 0)) = [];
    for k = 1:2:(filter_order)-1
        d = d + 1;
        polynomial = poly([P(k) P(k+1)]);
        C2 = 2/(R*polynomial(2));
        C5 = 1/(R^2*C2*polynomial(3));
        if lang == "spanish"
            fprintf('\n Par %d:\n Condensador 2: %e\n Condensador 5: %e\n',d,C2,C5)
        elseif lang == "english"
            fprintf('\n Par %d:\n Capacitor 2: %e\n Capacitor 5: %e\n',d,C2,C5)
        end
    end
    ugly_pole = -1*ugly_pole;
    syms C
    C = double(solve(ugly_pole == 1/(R*C)));
    syms R_fixed
    C_fixed = 100e-9;
    Rps = double(solve(ugly_pole == 1/(R_fixed*C_fixed)));
    syms R2
    R2 = double(solve(ugly_pole == (1/(R*C))*(1+(R2/R1))));
    if lang == "spanish"
        fprintf('\n Par %d:\n Resistencia R2: %e\n Condensador C: %e\n',d+1,R2,C)
        fprintf(' Resistencia Rps: %e\n',Rps)
    elseif lang == "english"
        fprintf('\n Par %d:\n Resistencia R2: %e\n Capacitor C: %e\n',d+1,R2,C)
        fprintf(' Resistor Rps: %e\n',Rps)
    end
end

%%Informative note
%This part of code is just for compile .m files in batch mode. 
%pause(15000)


