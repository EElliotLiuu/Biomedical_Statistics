# 生物医学统计基础
## Unit 1

1.生物医学数据类型：
- 数值型：离散/连续
- 类别型：有序分类/无序分类

2.生物医学统计常用术语
- 总体：研究对象的全体，用希腊字母表示其参数
- 样本：通常是总体的一部分，用英文字母表示统计量

3.推断性统计：用样本来推测总体
描述性统计：用全体的数字特征来描述总体

4.大数定理和中心极限定理
- 大数定理（LLN）：样本量足够大时，统计量接近于总体参数：
$$(m,s,r) \xrightarrow{n \to \infty } (\mu, \sigma, \rho)$$
- 中心极限定理（CLT）：无论总体是否符合正态分布，当样本量足够大时，样本均值近似服从正态分布。更特别的，若总体符合正态分布 $N(\mu, \sigma^2)$ ，则样本均值 $\bar{x} \sim N(\mu, \frac{\sigma^2}{n})$ 。


5.常用分布与函数

常用分布：

- 正态分布：$N(\mu, \sigma^2)$，Z分布为标准正态分布 $Z \sim N(1, 1)$

- 卡方分布：$X_i \stackrel{\text{i.i.d.}}{\sim} N(0,1),\quad \sum_{i=1}^{n}X_i^2 \sim \chi^2(df=n)$

- t分布：$Z \sim N(0, 1), \chi^2 \sim \chi^2(n), Z \perp \chi^2 \implies \frac{Z}{\sqrt{\chi^2/n}} \sim t(df=n)$

- f分布：$\chi_1^2 \sim \chi^2(m), \chi_2^2 \sim \chi^2(n), \chi_1^2 \perp \chi_2^2 \implies \frac{\chi_1^2/m}{\chi_2^2/n} \sim F(df_1=m, df_2=n)$


连续型概率分布通用函数
- **pdf（概率密度函数）**：描述连续型随机变量在取值x处的概率密集程度。连续随机变量单点概率恒为0，某区间内的概率由密度函数在该区间的积分得到。
- **cdf（累积分布函数）**：计算随机变量取值小于等于x的累积概率，即 $P(X \le x)$。
- **ppf（分位数函数）**：累积分布函数cdf的逆函数，输入累积概率q，返回对应的分位数值x（满足 $P(X \le x)=q$）。
- **sf（生存函数，又称残存函数）**：计算随机变量取值大于x的概率，即 $P(X > x) = 1 - \text{cdf}(x)$；在计算极小尾部概率时，数值精度显著高于直接用1减去cdf。
- **isf（逆生存函数）**：生存函数sf的逆函数，输入右侧尾部概率q，返回对应的上侧临界值x（满足 $P(X > x)=q$）。
- **rvs（随机数生成函数）**：从指定分布中抽取指定数量的随机样本。

离散型概率分布通用函数（以二项分布等计数分布为例）
离散分布同样包含cdf、ppf、sf、isf、rvs函数，语义与连续分布一致，仅输入取值为整数；同时用pmf替代连续分布的pdf：
- **pmf（概率质量函数）**：计算离散型随机变量恰好取某个整数值k的概率，即 $P(X = k)$。
- **cdf（累积分布函数）**：计算随机变量取值小于等于k的累积概率，即 $P(X \le k)$。
- **sf（生存函数）**：计算随机变量取值大于k的概率，即 $P(X > k) = 1 - \text{cdf}(k)$，处理极小概率时比直接做减法精度更高。
- **rvs（随机数生成函数）**：生成指定数量的随机样本，返回每个样本对应的离散取值（如二项分布的成功次数）。

6.数据可视化
- 看变化趋势：lineplot
- 看分布特征：boxplot，scatterplot，各种直方图。密度函数图
- 看关系：回归线，scatterplot+回归线，热图

## Unit 2

1.箱型图认读
- Median：中位数，50%分位数
- Q1：25%分位数
- Q3：75%分位数
- IQR=Q1-Q3
- $ Minimun=Q1-1.5 \times IQR $
- $ Maximun=Q3+1.5 \times IQR $

小于 Minimun 或大于 Maxinum 的为异常值

2.数据的集中程度

- 平均值 mean （对异常值敏感）
有n个观测值的样本$\{x_i\vert i=1,2,\ldots,n\}$定义
算术平均：$\overline{x}=\frac{x_1+x_2+\ldots+x_n}{n}$
几何平均（平均倍数变化率）：$\overline{x_g}=\sqrt[n]{x_1 x_2\ldots x_n}$
调和平均（平均比率）：$\overline{x_H}=\frac{n}{\frac{1}{x_1}+\frac{1}{x_2}+\ldots+\frac{1}{x_n}}$

- 众数 mode （用于离散型变量）
- 中位数 median （对异常值不敏感）

3.描述数据的分散性

极差 $Range=x_{(n)}-x_{(1)}$
变异系数 $cv=\frac{\sigma}{\mu}=\frac{\mathcal{s}}{\mathcal{m}}$
方差 $\sigma^2=\frac{\sum(x_i-\mu)^2}{N}$，$\mathcal{s}^2=\frac{\sum(x_i-\mathcal{m})^2}{n-1}$
标准差 $\sigma,\mathcal{s}$  

4.描述数据的分布

峰度 kurtosis
偏度 skewness
变异系数 coefficient of variation



**峰度**
Pearson's Kurtosis:
$$Kurtosis(x)=E[(\frac{x-\mu}{\sigma})^4]=\beta_2=\frac{\mu_4}{\sigma^4}$$
对于正态分布，上式恒等于3
Fisher's Kurtosis/Excess of Kurtosis = Pearson's Kurtosis - 3
FK=0时，为高斯分布

Fk>0时，为超高斯分布

FK<0时，为欠高斯分布

**偏度**
$$Skewness=\widetilde{\mu_3}=E[(\frac{X-\mu}{\sigma})^3]=\frac{\mu_3}{\sigma^3}$$ 

向左偏、偏度<0，向右偏、偏度>0（看尾巴）
<font color=red>注意，当偏度=0时仅能说明对称分布，不能说明正态性</font>

**变异系数**
$$CoV=\frac{\sigma}{\mu}\times 100\%=\frac{\mathcal{s}}{\mathcal{m}}\times 100\%$$
用途：  

1. 比较不同数据集的离散程度
2. 无量纲
缺点：当均值接近0时不好用  

## Unit 3
1.统计推断三大方法
- 置信区间
- 零假设显著性检验
- 贝叶斯统计（Not required）


2.置信区间

对于单个样本的均值，常用的枢轴量有：
- **总体方差 $\sigma^2$ 已知**：$\bar{x} \sim N(\mu, \frac{\sigma^2}{n})$, 标准化后即 $Z = \frac{\bar{x}-\mu}{\sigma/\sqrt{n}} \sim N(0, 1)$。
- **总体方差 $\sigma^2$ 未知**：采用样本方差 $s^2$，则枢轴量为 $T = \frac{\bar{x}-\mu}{s/\sqrt{n}} \sim t(n-1)$。

配对样本 (Paired Samples)
来自同一对象的两次测量（如药前/药后）。通过构造差值 $d_i=x_{1i}-x_{2i}$ 转化为单样本情况，计算 $\bar{d}$ 和 $s_d$，自由度为 $n-1$。

样本均值差异的置信区间

总体标准差已知

当 $\sigma_1$ 和 $\sigma_2$ 已知时，我们可以知道两个样本均值 $x_1 \sim N(\mu_1, \frac{\sigma_1^2}{n_1})$，$x_2 \sim N(\mu_2, \frac{\sigma_2^2}{n_2})$，所以其样本均值之差:

 $$
 \bar{x_1}-\bar{x_2} \sim N(\mu_1-\mu_2, \frac{\sigma_1^2}{n_1}+\frac{\sigma_2^2}{n_2})
 $$

标准误 $SE=\sqrt{\frac{\sigma_1^2}{n_1}+\frac{\sigma_2^2}{n_2}}$，所以置信区间为：

$$
CI=(\bar{x}_1-\bar{x}_2)\pm Z_{\alpha/2} \cdot \sqrt{\frac{\sigma_1^2}{n_1}+\frac{\sigma_2^2}{n_2}}
$$

总体标准差未知但相等

当 $\sigma_1=\sigma_2=\sigma$但未知时，我们采用合并方差法：$s_p^2=\frac{(n_1-1)s_1^2+(n_2-1)s_2^2}{n_1+n_2-2}$，此时我们构造T统计量：

$$
T = \frac{(\bar{X}_1 - \bar{X}_2) - (\mu_1 - \mu_2)}{s_p \cdot \sqrt{\frac{1}{n_1} + \frac{1}{n_2}}} \sim t(n_1+n_2-2)
$$

标准误 $SE=s_p \cdot \sqrt{\frac{1}{n_1} + \frac{1}{n_2}}$，所以置信区间为：

$$
CI=(\bar{x}_1-\bar{x}_2)\pm t_{\alpha/2,n_1+n_2-2} \cdot s_p \sqrt{\frac{1}{n_1}+\frac{1}{n_2}}
$$

总体标准差未知且不等

这是最普通的情况，当 $\sigma_1 \neq \sigma_2$ 时，我们不再合并方差，而是直接使用各自的样本方差。此时我们构造T统计量：

$$
T = \frac{(\bar{X}_1 - \bar{X}_2) - (\mu_1 - \mu_2)}{\sqrt{\frac{s_1^2}{n_1} + \frac{s_2^2}{n_2}}} \sim t(df)
$$

标准误 $SE=\sqrt{\frac{s_1^2}{n_1} + \frac{s_2^2}{n_2}}$，而此时的自由度 $df \neq n_1+n_2-2$，而需要用Welch-Satterthwaite 自由度公式来计算：

$$
df=\frac{\left( \frac{s_1^2}{n_1}+\frac{s_2^2}{n_2} \right) ^2}{\frac{1}{n_1-1}\left(\frac{s_1^2}{n_1}\right)^2
+\frac{1}{n_2-1}\left(\frac{s_2^2}{n_2}\right)^2}
$$

最终求得置信区间为：

$$
CI=(\bar{x}_1-\bar{x}_2)\pm t_{\alpha/2,df} \cdot \sqrt{\frac{s_1^2}{n_1}+\frac{s_2^2}{n_2}}
$$