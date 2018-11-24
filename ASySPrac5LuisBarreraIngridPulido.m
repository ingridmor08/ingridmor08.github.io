%% Práctica 5: Series de Fourier en tiempo continuo
%
% 
%% Integrantes
%
% Barrera Bautista Luis Franciso
%
% Pulido Morales Ingrid 
%% Objetivos
%
% * Realizar gráficas de series de Fourier exponenciales y trigonómetricas
% en tiempo coontinuo
% * Manipulación de instrucciones en MATLAB
% * Calculo númerico de los coeficientes de Fourier
%% Introducción
%
% Aproximación Numérica de $D_{n}$
%
% Podemos calcular $D_{n}$ usando el DTF (Transformada Discreta de
% Fourier), la cual usa las muestras de una señal periódica $x(t)$ en un
% periodo. El intervalo de muestreo es T segundos. Por lo tanto, hay
% $N_{0}=\frac{T_{0}}{T}$ numero de muestras en un periodo $T_{0}$. Para
% encontrar una relación entre $D_{n}$ y las muestras de $x(t)$,
% consideramos
% 
% $D_{n}=\frac{1}{T_{0}}\int_{T_{0}}x(t)e^{-jnw_{0}t}dt$
%
% $=\lim_{T\rightarrow 0} \frac{1}{N_{0}T}\sum_{k=0}^{N_{0}-1}x(kT)e^{-jnw_{0}kT}T$
%
% $=\lim_{T\rightarrow 0} \frac{1}{N_{0}}\sum_{k=0}^{N_{0}-1}x(kT)e^{-jn\Omega k} \ \ \ \ \ \ \ \ \ \ \ \ (1)$
%
% Donde $x(kT)$ es la kth muestra de $x(t)$ y
%
% $N_{0}=\frac{T_{0}}{T} \ \ \ \ \ \ \ \ \ \Omega_{0}=\omega_{0}T= \frac{2\pi}{N_{0}}\ \ \ \ \ \ \  \ \ \ \ \ \ \ (2)$ 
%
% Note que es imposible hacer $lim_{T\rightarrow 0}$ en la Ec.(1). Podemos hacer $T$ pequeña, pero no hacerla cero, lo que causará que los datos se incrementen sin un límite.
% Entonces, ignoramos el límite en $T$ en la Ec. (1) con el argumento de
% que $T$ es razonablemente pequeño. Por lo tanto, podemos expresar la Ec. (1) como: 
%
% $D_{n}=\frac{1}{N_{0}}\sum_{k=0}^{N_{0}-1} x(kT)e^{-j\pi\Omega_{0}k}\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (3)$
%  
% Ahora, de la Ec.(2), $\Omega_{0}N_{0}=2\pi$. Por lo tanto,
% $e^{j\Omega_{0}(k+N_{0})}$ $=e^{j\Omega_{0}k}$ y de la Ec. (3)
%
% $D_{n+N_{0}}=D_{n}$  
%
% La propiedad de periodicidad $D_{n+N_{0}}=D_{n}$ significa que más allá
% de $n=\frac{N_{0}}{2}$, los coeficientes representan los valores para las
% $n$ negativas. Por ejemplo $N_{0}=32,
% D_{17}=D_{-15},D_{18}=D_{-14},...,D_{31}=D_{-1}$. El ciclo se repite de
% nuevo en n=32.
%% Desarrollo
%% Ejemplo 6.1
% $f(t)= e^{-\frac{t}{2}}$
% 
% Se calculó la Serie de Fourier trigonométrica, popr lo tanto sus
% coeficientes son
%
% $a_{n}=\frac{2}{\pi}\int_{0}^{\pi}e^{-\frac{t}{2}}cos(2nt)dt= 0.504(\frac{2}{1+16n^{2}})$
%
% $b_{n}=\frac{2}{\pi}\int_{0}^{\pi}e^{-\frac{t}{2}}sin(2nt)dt=0.504(\frac{8n}{1+16n^{2}})$
%
% Graficando la serie
%
% * Para 4 armónicos 
%
% <<ej1_1.png>>
%
% * Para 15 armónicos
%
% <<ej1_2.png>>
% 
%% Ejemplo 6.2
%
% <<fig1.PNG>>
%
% Se calculó la Serie de Fourier Exponencial 
%  
% $D_{n}=\frac{1}{2}\int_{-0.5}^{1.5}x(t)e^{-jn \pi t} = -\frac{12jsin(\frac{n \pi}{2})}{n^{2}\pi^{2}}$
%
% Graficando la serie
%
% * Para 4 armónicos
%
% <<ej2_1.png>>
%
% * Para 15 armónicos
%
% <<ej2_2.png>>
%
%% Ejemplo 6.4
%
% Se calculó la Serie de Fourier Exponencial para la siguiente función
%
% <<fig2.PNG>>
%
% Calculando su $D_{n}$
%
% $D_{n}=\frac{1}{2\pi}\int_{0}^{\pi}e^{-jnt}dt= \frac{sin(\frac{\pi n}{2})}{\pi n}$
%
% Graficando la Serie
%
% * Para 4 armónicos
%
% <<ej3_1.png>>
%
% * Para 15 armónicos
%
% <<ej3_2.png>>
%% Ejercicio 6.5
%
% $x(t)=\left | sin(t) \right |$ 
%
% $D_{n}=\frac{1}{\pi}\int_{0}^{\pi}\left | sin(t) \right | e^{-jn\pi t}dt = \frac{2}{\pi - 4\pi n^{2}}$
%
% Graficando la serie
%
% * Prara 4 armónicos
%
% <<eje4_1.png>>
%
% * Para 15 armónicos
%
% <<eje4_2.png>>
%% Ejemplo 6.7
%
% $x(t)=\delta_{T_{0}}(t)$
%
% Calculando el $D_{n}$ con $T_{0}= 3$
%
% $D_{n}= \frac{1}{3} \int_{-\frac{3}{2}}^{\frac{3}{2}} \delta (t)e^{-jn\omega_{0}t}dt$
%
% Graficando la serie
%
% * Para 4 armónicos
%
% <<eje5_1.png>>
%
% <<eje5_2.png>>
%
%% Punto 6
%
% Se elaboró un código similar al COMPUTER EXAMPLE C6.2 que se muestra a
% continuación
% 
%   t0=-0.5;
%   tf=1.5;
%   t=linspace (t0, tf,1000);
%   x = @(t) 6.*t.*(t>=-0.5 & t<0.5)+6.*(1-t).*(t>=0.5 & t<1.5);
%   t = linspace (-2*pi, 2*pi,1000);
%   y=x(t);
%   sumterms = zeros(16, length(t));
%   gei Ilinspace();
%   sumterms(1,:) = 1/2;
%   for n = 1:size(sumterms,1)-1;
%   sumterms(n+1,:) = 24*sin(n*pi/2)*sin(pi*n*t)/(n^2*pi^2);
%   end
%   x_N = cumsum (sumterms); 
%   figure(1); 
%   clf; 
%   ind = 0;
%   for N = [0,1:2:size(sumterms, 1)-1]
%   ind = ind+1; 
%   subplot (3,3,ind);
%   plot (t,x_N(N+1,:),t,y,'k--'); 
%   axis ([-1 2 -4 4]);
%   xlabel('t');
%   aux=strcat('x_{',num2str(N),'}');
%   ylabel(aux);
%   end
%
% <<eje6.png>>
%
%% Punto 7
%
% Para el ejemplo 6.1 se implementaron el algoritmo de trapecio compuesto y
% el código COMPUTER EX C6.4 para calcular desde $D_{0}...D_{4}$ con n=15
%
% Realizando una tabla de comparación
%
% <<fig3.PNG>>
%
%% Bibliografía
%  B. P. Lathi. (2005). Linear Systems and Signals. Oxford, EE. UU.: Oxford University Press.
%  Martinez, M. Rafael. (2018). Series de Fourier para señales continuas. 2018, de Mate y así Sitio web: https://www.youtube.com/user/rafa5131
%%