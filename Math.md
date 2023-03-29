# 分子模拟中用到的数学公式

## LJ势与牛顿力学之间的关系

​		分子动力学中, 最基本最简单的势函数就是LJ势, 而LJ势也有很多衍生种, 这里我们介绍的是最基本的截断LJ势, 具体的数学公式表达如下
$$
U_i\left(r_{i j}\right)=\left\{\begin{array}{ll}
4 \epsilon\left[\left(\frac{\sigma}{r_{i j}}\right)^{12}-\left(\frac{\sigma}{r_{ij}}\right)^{6}\right] & r_{i j}<r_{c} \\
0 & r_{i j} \geq r_{c}
\end{array}\right. \tag{1.1}
$$

​		其中$r_c$为LJ势的截断半径, 而$r_{ij}$为第i个原子和第j个原子之间的距离, $U_i$为第i个原子之间的势能, $\epsilon$和$\sigma$是LJ势中计算势能要用到的参数

​		而在本例中, 我们将通过忽略有吸引力的尾部并将上式改为
$$
U_i\left(r_{i j}\right)=\left\{\begin{array}{ll}
4 \epsilon\left[\left(\frac{\sigma}{r_{i j}}\right)^{12}-\left(\frac{\sigma}{r_{i j}}\right)^{6}\right]+\epsilon & r_{i j}<r_{c}=2^{1 / 6} \sigma \\
0 & r_{i j} \geq r_{c}
\end{array}\right. \tag{1.2}
$$


​		使用该电位构建的模型流体只不过是一组既柔软（尽管柔软度有限）又光滑的碰撞球。所有将系统连接在一起的是容纳原子（或球）的容器。虽然这种高度简化的模型可以定量表示的系统种类有限，通常是低密度气体，但它与更为精确的模型有很多共同点，并且在计算简单性方面具有明显的优势。

​		那么我们同理可以获得计算原子之间受力的公式通过经典力学
$$
\vec{f}=-\nabla u(\vec{r}) \tag{1.3}
$$
​		将其展开, 可以获得第i个和第j个原子之间的作用力为
$$
f_{i j}=\left(\frac{48 \epsilon}{\sigma^{2}}\right)\left[\left(\frac{\sigma}{r_{i j}}\right)^{14}-\frac{1}{2}\left(\frac{\sigma}{r_{i j}}\right)^{8}\right]r_{ij} \boldsymbol  ~~~~~~~~~~~~~~~ (r_{ij}<r_{c}) \tag{1.4}
$$
​		当然(4)式的成立前提条件一定不能忘记是小于截断半径, 再利用牛顿第二定律, 有
$$
m \vec{\ddot{\boldsymbol{r}_{i}}}=\vec{{f}_{i}}=\sum_{\substack{j=1 \\(j \neq i)}}^{N_{m}} \vec{f_{i j}} \tag{1.5}
$$
​		$m$是原子质量, $\vec{\ddot{\boldsymbol{r}}_{i}}$ 为第i个原子的加速度, $N_{m}$为系统中所有原子的数目, 同时注意式(5)是矢量计算, 根据牛顿第三定律, $\vec{f_{ij}} = \vec{-f_{ij}}$ ,因此每个原子对仅需要计算一次即可, 那么正常来说, 需要计算$\frac{N_{m}(N_m - 1)}{2}$次的循环即可完成所有原子的受力计算, 而上述算法也被称为verlet算法.

​		通过适当的转换单位, 可以简化LJ势的表达式, 减小计算量, 并提供数值计算的精度, 这里提供一组LJ势下的单位转换:

​													长度: $r \rightarrow r\sigma $ 

​													能量: $e \rightarrow e\epsilon $

​													时间: $t \rightarrow t \sqrt{m\sigma^2\epsilon} $

​		则运动方程中, 加速度的表达式可以改写为:
$$
\bold {\vec a_i} = 48 \sum_{j(\neq i)}({\bold r_{ij}}^{-12} - \bold{r_{ij}}^{-6}) \bold r_{ij} \tag{4.1}
$$
 同时, 一整个系统的势能和动能可以按照以下公式进行计算:
$$
E_{K}=\frac{1}{2 N_{m}} \sum_{i=1}^{N_{m}} {v}_{i}^{2} \tag{1.6}
$$

$$
E_{U} = \frac{4}{N_m} \sum_{1\leq i \lt j \leq N_m} ({r_{ij}^{-12}} - {r_{ij}^{-6}}) \tag{1.7}
$$ {E}

​		$E_k$和$E_U$分别代表整个系统的动能和势能, 而温度的计算形势与上述的动能计算性能很像, 因为我们知道温度和原子的运动剧烈程度有关, 原子运动越剧烈, 温度越高
$$
T = \frac{1}{dN_m} \sum_i{{v_i}^2} \tag{1.8} = 2 \frac{E_k}{d}
$$
​		其中d为维数(dimension), 取值为2或者3, 根据自己的仿真情况而定, 其余参数与上述公式相同. 这里我们可以发现对于二维的情况($d= 2$), $T$ 正好等于 $E_k$.

## Taylor展开

​		在聊到蛙跳积分法之前, 我们先来聊一聊泰勒展开, 对于一个函数$f(t)$来说, 想要求取$f(t+ \Delta t)$的值, 根据高等数学的定义, 有如下的等式成立:
$$
f(t + \Delta t) = \int_{t}^{t + \Delta t} {f(x)}' dx + f(t) \tag{1}
$$
​		现在, 我们将目光放到积分项上面, 根据分部积分公式有:
$$
\int_{t}^{t + \Delta t} {f'(x)} dx = ((x+C){f(x)}')|_{t}^{t+\Delta t} - \int_{t}^{t + \Delta t}(x+C){f''(x)}dx \tag{2}
$$
​		在上述转换过程中:

​									               $u = f'(x) ~~~~~~~~~~~~~~~~~~~~~~~~~~ dv = dx$

​											$du = f''(x)dx ~~~~~~~~~~~~~~~~~~~~~~~~~~v= x + C$ 

​		这里我们可以设$C = -(t + \Delta t)$, 那么(2)式就可以改写为如下形势:
$$
\int_{t}^{t + \Delta t} {f'x)} dx = ((x-(t + \Delta t)){f(x)}')|_{t}^{t+\Delta t} - \int_{t}^{t + \Delta t}(x-(t + \Delta t)){f'x)}dx \tag{2}
$$
​		这是因为在$x = t + \Delta t $的时候, $ ((x-(t + \Delta t)){f(x)}') =0$, 那么化简上式之后, 我们可以获得:
$$
\int_{t}^{t + \Delta t} {f(x)}' dx = \Delta tf'(x) -\int_{t}^{t + \Delta t}(x-(t + \Delta t)){f(x)}''dx \tag{3}
$$
​		同理, 对于(3)中的积分项, 我们可以继续使用分部积分法进行展开, 这里不再赘述, 直接给出最终形式
$$
f(t + \Delta t) = f(t) + \Delta tf'(x) + \Delta t^2 \frac{f''(x)}{2! } + \Delta t^3 \frac{f'''(x)}{3! } + \dots \tag{4}
$$

$$
f(t - \Delta t) = f(t) - \Delta tf'(x) + \Delta t^2 \frac{f''(x)}{2! } - \Delta t^3 \frac{f'''(x)}{3! } + \dots \tag{5}
$$

​		将(4)式和(5)式合并之后可得:
$$
f(t + \Delta t) = 2f(t) - f(t - \Delta t) \tag{6} + \Delta t^2 f''(x) 
$$
​		根据Taylor展开的定义, 我们不难证明, 此时的$f(t + \Delta t)$具有$o(\Delta t^4)$的精度. 至此, leapfrog法的理论基础已经叙述完毕.

## leapfrog法积分(蛙跳积分法)

​		首先我们定义$h = \Delta t$为一个最小的时间积分步长, 根据泰勒展开, 可以按以下公式更新原子坐标和速度
$$
v_{ix}(t + \frac{h}{2}) = v_{ix}(t - \frac{h}{2}) + ha_{ix}(t) \tag{2.1}
$$

$$
r_{ix}(t + h) = r_{ix}(t) + hv_{ix}(t + \frac{h}{2}) \tag{2.2}
$$

​		上述算法可以看到位置和速度不是在同一个时刻更新, 如果想要在同一个时刻更新速度, 则需要:
$$
v_{ix}(t) = v_{ix}(t - \frac{h}{2}) + {\frac{h}{2} }a_{ix}(t) \tag{2.3}
$$
​		对于y方向上的速度与坐标更新处理也是同理.

​		蛙跳法可以用另一种代数等价的方式重新公式化，使坐标和速度能够在同一时间点进行评估，从而避免了（2.3）中的速度调整。为此，计算分为两部分：在计算加速度值之前，使用旧的加速度值将速度更新半个时间步长，
$$
v_{i x}(t+h / 2)=v_{i x}(t)+(h / 2) a_{i x}(t) \tag{2.4}
$$

$$
r_{i x}(t+h)=r_{i x}(t)+h v_{i x}(t+h / 2) \tag{2.5}
$$

​		然后使用中间速度值将坐标更新一个完整的时间步长
$$
v_{i x}(t+h)=v_{i x}(t+h / 2)+(h / 2) a_{i x}(t+h)\tag{2.6}
$$
​		基于上述积分法计算出来位置$\vec{r}$具有$o(h^4)$的精度, 而速度$\vec{v}$具有$o(h^2)$的精度.



## 基本物理量的测量

​		在本案例研究中，能量和压力是唯一测量的属性。压力根据维里表达式确定:
$$
PV = N_mT + \frac{1}{d} \left\langle\sum_{i=1}^{N_m} \boldsymbol{r_i} \cdot \boldsymbol{f_i} \right\rangle \tag{3.1}
$$
​		在二维计算中, 体积$V$可以用区域的面积来代替, 对于势函数来说, 可以将其写作如下形势表明其遍历了每一个原子对:
$$
PV = N_mT + \frac{1}{d} \left\langle\sum_{i \lt j} \boldsymbol{r_i} \cdot \boldsymbol{f_i} \right\rangle \tag{3.2}
$$
​		将LJ势的表达式带入其中可得:
$$
PV = N_mT + \frac{1}{d} \left\langle\sum_i v_i^2 + 48 \sum_{i \lt j} (r_{ij}^{-12} - \frac{1}{2}r_{ij}^{-6})  \right\rangle \tag{3.1}
$$
​		同时总能量$E=E_K + E_U$理论上在仿真过程中应该总是守恒的.





