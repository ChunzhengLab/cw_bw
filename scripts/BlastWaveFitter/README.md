# Blast‐Wave Fitter 公式集
### Chunzheng Wang

## 1 流速剖面（含二阶调制）
$$
\rho(r,\phi)\;=\;r^{\,n}\,\Bigl[\rho_0+\rho_2\cos\bigl(2\phi_b\bigr)\bigr],
\qquad r\in[0,1].
$$
其中 $n$ 为可调的径向幂指数。

- 若使用椭圆冻结面，则
$$
x = R_x\,r\cos\phi,\quad
y = R_y\,r\sin\phi,
$$

面积元为 $dA = R_x R_y\, r\,dr\,d\phi$。

## 1.1 坐标与角度变换
坐标映射公式已在上文第 1 节给出，故此处省略复述。

- 流速方向角 $\phi_b$ 与位置角 $\phi$ 的关系  
  $$
    \tan\phi_b
    =\frac{R_x^2}{R_y^2}\tan\phi
  $$

下面所有含 $\phi_b$ 的余弦项与权重均使用这一角度。

## 1.2 平均径向流速与快速度的互换

- 外缘快速度 $\rho_0$ 与流速 $\beta_T$ 之间的转换  
  $$
    \rho_0 = \operatorname{atanh}(\beta_T),\qquad
    \beta_T = \tanh(\rho_0).
  $$

## 1.3 二阶调制参数的数学意义

- 径向快速度剖面  
  $$
    \rho(r,\phi)
    = r^{\,n}\,\bigl[\rho_0+\rho_2\cos(2\phi_b)\bigr].
  $$

## 2. 单粒子谱 (Spectrum)

本模型假设在每个流体元的局部共动系中，粒子的能量分布服从**相对论 Maxwell–Jüttner 分布**：
$$
  f_0(E^*) \;\propto\; p^*\,E^*\,\exp\!\Bigl(-\frac{E^*}{T_{\rm kin}}\Bigr),
  \quad E^* = p^\mu u_\mu = m_T\cosh\rho - p_T\sinh\rho\cos(\phi_p-\phi_b).
$$
在对该分布在角度和径向坐标上做积分并归一化后，即得到下面所示的包含修正贝塞尔函数 $K_1$ 和 $I_0$ 的谱表达式。

$$
\frac{1}{2\pi p_T}\frac{d^2N}{dp_T\,dy}
=\;N_0 R_x R_y \!\int_{0}^{1}\!r\,dr\!\int_{0}^{2\pi}\!d\phi\;
m_T\,K_1\!\Bigl(\frac{m_T\cosh\rho}{T_{\rm kin}}\Bigr)\,I_0\!\Bigl(\frac{p_T\sinh\rho}{T_{\rm kin}}\Bigr),
$$
其中 $ m_T = \sqrt{p_T^2 + m^2} $, $ K_1(z) $ 为第二类修正贝塞尔函数，$I_\nu(z)\ (\nu=0,2)$ 为第一类修正贝塞尔函数，分别对应谱和 $v_2$ 中的 $I_0,I_2$.


## 3. 椭圆流系数 $v_2$

- 约定质量横截面 $m_T$  
  $$
    m_T = \sqrt{p_T^2 + m^2}.
  $$

- 定义广义积分  
  $$
    I_\nu(p_T)
    = R_x R_y \!\int_0^1 r\,dr \int_0^{2\pi} d\phi\;
      m_T\,K_1\!\Bigl(\tfrac{m_T\cosh\rho}{T_{\rm kin}}\Bigr)\,
      I_\nu\!\Bigl(\tfrac{p_T\sinh\rho}{T_{\rm kin}}\Bigr),
    \quad \nu=0,2.
  $$
  （其中对 $\nu=2$ 还需乘 $\cos2\phi_b$）

- 谱与流统一表达  
  - 谱：$\nu=0$ 对应 $I_0$。  
  - $v_2(p_T)$ 分子取 $\nu=2$ 并乘 $\cos2\phi_b$，分母取 $I_0$，即  
    $$
      v_2(p_T)
      = \frac{\displaystyle
        \int r\,dr\,d\phi\; m_T K_1\!\bigl(\tfrac{m_T\cosh\rho}{T_{\rm kin}}\bigr)\,
        I_2\!\bigl(\tfrac{p_T\sinh\rho}{T_{\rm kin}}\bigr)\cos2\phi_b}
      {\displaystyle
        \int r\,dr\,d\phi\; m_T K_1\!\bigl(\tfrac{m_T\cosh\rho}{T_{\rm kin}}\bigr)\,
        I_0\!\bigl(\tfrac{p_T\sinh\rho}{T_{\rm kin}}\bigr)}.
    $$
$$
v_2(p_T)
=\frac{\displaystyle
R_x R_y\!\int_{0}^{1}\!r\,dr\!\int_{0}^{2\pi}\!d\phi\;
m_T\,K_1\!\Bigl(\frac{m_T\cosh\rho}{T_{\rm kin}}\Bigr)\,I_2\!\Bigl(\frac{p_T\sinh\rho}{T_{\rm kin}}\Bigr)\cos\bigl(2\phi_b\bigr)
}
{\displaystyle
R_x R_y\!\int_{0}^{1}\!r\,dr\!\int_{0}^{2\pi}\!d\phi\;
m_T\,K_1\!\Bigl(\frac{m_T\cosh\rho}{T_{\rm kin}}\Bigr)\,I_0\!\Bigl(\frac{p_T\sinh\rho}{T_{\rm kin}}\Bigr)
}.
$$

## 4. 参数定义

- $\beta_T$: 最大径向流速（无量纲，$0<\beta_T<1$）。
- $T_{\rm kin}$: 运动学冻结温度（GeV）。
- $n$: 流速剖面指数。
- $\rho_2$: 二阶各向异性振幅。
- $N_0$: 谱线归一化常数。
- $R_x, R_y$: 椭圆冻结面半长轴，通常固定$R_y=10\,\mathrm{fm}$，仅拟合$R_x$。
