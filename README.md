# SPE10


## Задача 2-фазной фильтрации: метод Ньютона}

### Постановка задачи:


$
\label{trivial}
\left\{\begin{matrix}
\phi \frac{\partial S}{\partial t} - div(S_{k}\bigtriangledown p) = g_{0}=0
\\
\phi \frac{\partial (1-S)}{\partial t} - div((1-S)k\bigtriangledown p) = g_{\omega }=0
\end{matrix}\right.
$

 Рассматривается скважина радиуса r. T - полное время интегрирования по времени, $\bigtriangleup t$ - шаг по времени, $t_0 = 0$, m - количество шагов, $S\in \left [ 0,1 \right ]$ - насыщение нефти, $S_0 = 0.75$, p - давление, $p_0 = 100$, $\phi = const$- пористость.
 
 Схема решения:
 
 0) Нулевое Ньютоновское приближение: $p^{0}=p_{0}$, $S^{0}=S_{0}$
 
 1) Изначально принимается $k=1$. По формуле невязки считаем $R_{0}^{k}$ и $R_{\omega}^{k}$:
$
\label{trivial}
R_{0}=\phi \frac{\partial S}{\partial t}-div(Sk\bigtriangledown p)-q_{0}=\phi \frac{S_{i}-S_{i}^{0}}{\bigtriangleup t}+\sum_{j}T_{ij}\begin{pmatrix}S_{i}, p_{i}> p_{j}
\\ S_{j}, p_{i}<  p_{j}
\end{pmatrix}-k(S_{i})WI(p_{bh}-p_{i})=0,
$

где $k(S_i) = 1$, если нагнетающая скважина. Если добывающая, то $k(S_i) = S_i$

2) Проверка условия:

$
\label{trivial}
\left \|R_{0}^{k}  \right \|+\left \|R_{\omega }^{k}  \right \|< \varepsilon 
$

если условие выполнено, то $p_{0}=p^{k}$ и $S_{0}=S^{k}$. Если нет, переходим к 3 пункту.

3) Вычисляем Якобиан:

$
\label{trivial}
J_{0}=\left (\frac{\partial R_{0} }{\partial p^{T}} \frac{\partial R_{0} }{\partial S^{T}} \right )
$

$
\label{trivial}
J_0=\left\{\begin{matrix}\frac{\partial R_{0i} }{\partial p_{i}} =\sum T_{ij}\begin{pmatrix}S_{i}, p_{i}> p_{j}
\\ S_{j}, p_{i}<  p_{j}
\end{pmatrix}-k(S_{i})WI
\\ \frac{\partial R_{0i}}{\partial p_{j}}=-T_{ij}\begin{pmatrix}S_{i}, p_{i}> p_{j}
\\ S_{j}, p_{i}<  p_{j}
\end{pmatrix}
\\ \frac{\partial R_{0i}}{\partial S_{i}}=\frac{\phi }{\bigtriangleup t}+\sum T_{ij}(p_i - p_j)\begin{pmatrix}1, p_{i}> p_{j}
\\ 0, p_{i}<  p_{j}
\end{pmatrix}+WI(p_i-p_{bh})\frac{\partial k(S_i)}{\partial S_i}
\\ \frac{\partial R_{0i}}{\partial S_j}=T_{ij}(p_i-p_j)\begin{pmatrix}1, p_{i}> p_{j}
\\ 0, p_{i}<  p_{j}
\end{pmatrix}
\end{matrix}\right.
$


Используется 5-точечная схема конечных разностей для Якобиана:
 
$
\label{trivial}
\tau_1u_{i-1,j}-\tau_1u_{i,j}+\tau_2u_{i,j+1}-\tau_2u_{i,j}+\tau_3u_{i+1,j}-\tau_3u_{i,j}+\tau_4u_{i,j}
-\tau_4u_{i,j-1}=r(x)=0
$

<!-- \begin{figure}[h!] -->
  <!-- \begin{center} -->
  <!-- \scalebox{0.5}{ -->
<!-- \includegraphics{5.png} -->
 <!-- }  -->
 <!-- \end {center} -->
<!-- \end {figure} -->

4) Решаем систему:

$
\label{trivial}
J_0\begin{bmatrix}\bigtriangleup p
\\ 
\bigtriangleup S
\end{bmatrix} = \begin{bmatrix}-R_{0}^{k}
\\ 
-R_{\omega }^{k}
\end{bmatrix}
$

Если невязка сходится:
$
\label{trivial}
\begin{Vmatrix}-R_0
\\ 
-R_\omega 
\end{Vmatrix}=0, 
$

то задача решена: $p_{0}=p^{k}$ и $S_{0}=S^{k}$.

5) Следующий шаг Ньютона:

$
\begin{matrix}
p^{k+1}=p^{k}+\alpha\bigtriangleup p\\ 
S^{k+1}=S^{k}+\alpha\bigtriangleup S,
\end{matrix}
$

где $\alpha$ - подбирается так, чтобы значения давления и насыщенности не выходили за область допустимых значений. 

6) Шаг увеличивается: $k=k+1$. $m=m+1 \Rightarrow t = \bigtriangleup t*m$, когда $\bigtriangleup t*m\geq T$ алгоритм заканчивается на последней итерации, если условие не выполнено, то переходим к пункту 1

<p align="center">
  <img src="data/Example.png">
</p>