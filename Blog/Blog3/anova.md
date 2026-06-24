# 方差分析 ANOVA
基于零假设检验，我们引入 analysis of varience(ANOVA)。

不妨想象这样的一个现实问题：某个单次检验的假阳性率为 $\alpha$ ，对50个样本进行检验，则 $n$ 次独立检验后，出现假阳性的概率为 $FPR=1-(1-\alpha)^n$。以我们在NHST中常用的显著性水平 $\alpha=0.05$ 为例，若检测100次， 假阳性率将达到99.4%，非常惊人！

所以，如果我们在对多个样本进行两两配对NHST时，出现第一类错误的概率是累积的。为了避免以上这种情况，在检验多组均值是否存在差异时，我们先做整体检验，也就是ANOVA。

## ANOVA 的原理
### 1.lead-in
首先我们来看下面三个样本：三个样本的总体都符合正态分布，方差 $\sigma^2$ 均为1，不过均值 $\mu$ 分别为0、1、2，每个样本都是从总体随机采样100个数据得到的：

<div align="center">
  <img src="anova1.png" width="500">
</div>

可以看到，这三组数据之间并不好区分。如果我们把三个总体的方差设置小一点呢，比如让 $\sigma^2=0.1$，那么得到的分布图如下：

<div align="center">
  <img src="anova2.png" width="500">
</div>

可以看到，三组数据之间已经可以区分得很清楚了。哪还有没有其他的方法呢，我们让三个总体均值 $\mu$ 的区别变大一点，分别为0、4、8，那么得到以下的分布图：

<div align="center">
  <img src="anova3.png" width="500">
</div>

可以看到也比最开始的情况好区分。

### 2.ANOVA的数学模型
从上面的例子中可以看到，不同数据之间和组内的差异都对我们区分不同总体有影响。组间差距越大、组内差距越小，我们越能区分不同的总体。

抽象为数学模型：假设一共有 $k$ 组数据，分别属于不同的总体，每组数据都有 $n_i, i=1,2,\cdots,k$ 个观测值，则总共有 $ \sum_{i=1}^{k} n_i=N $ 个观测值，我们设 第 $i$ 组的第 $j$ 个观测值为 $x_{ij}$。

- 总平方和： $$ SST=\sum_{i=1}^{k}\sum_{j=1}{n_i}(x_{ij}-\bar{x}_{\cdot\cdot})^2$$

- 组间平方和：$$ SSB=\sum_{i=1}^{k}n_i(\bar{x}_{i\cdot}-\bar{x}_{\cdot\cdot})^2 $$

- 组内平方和：$$ SSW(SSE)=\sum_{i=1}^{k}\sum_{j=1}{n_i}(x_{ij}-\bar{x}_{i\cdot})^2 $$

- 根据数学公式推导，可知 $SST=SSB+SSW(SSE)$

- 组间方差 $$ MSB=\frac{SSB}{df_B}=\frac{SSB}{k-1} $$

- 组内方差 $$ MSW=\frac{SSW}{df_W}=\frac{SSW}{N-k} $$

关于自由度消失的问题，推导过程可以参考上一个Blog，组间方差只有一个总的平均值的约束，所以自由度减1；而组内方差有 $k$ 个组内平均值的约束，所以自由度减 $k$。

我们定义F统计量：$ F=\frac{MSB}{MSW} \sim F(df_B, df_W) $

效应量 $ \eta^2=\frac{SSB}{SST} $，表示SST中有多少比例是可由SSB解释的

可以看到，样本之间的差异越大，F统计量的值也越大。

### 3.ANOVA的类型
如果实验中同时考察了n个因素，我们称为 n-way ANOVA，n因素方差分析，因素一般可以分为组间因素和组内因素两种。
同时，如果我们对同一组对象同不同的试验，涉及重复测量，则称 RM n-way ANOVA，比如对一组病人，我们分析他们用了不同药物前后的变化，就涉及重复测量。
我们主要讨论 1-way ANOVA和 2-way ANOVA 两种情况。

## 1-way ANOVA
### 1.Python函数
首先，1-way ANOVA的零假设为 $H_0: \mu_1=\mu_2=\mu_3$，备选假设为 $ H_1: \mu_1=\mu_2=\mu_3 $ is not true，也就是三个均值至少有两个不相等。

然后是我们的anova函数，是在`pingouin`包里，一般在头文件

```py
import pingouin as pg
```

这里简述一下Python函数的用法，假如我们要对`data`数据进行 1-way ANOVA，就要使用函数：
```py
pingouin.anova(
  data=data, 
  dv='x', 
  between="Group", 
  ss_type=2, 
  detailed=True, 
  effsize='np2')
```
- `dv` 是因变量列名，比如在研究不同药物的作用时，dv这列就应该是作用效果
- `between` 是分组变量列名，单因素用字符串，多因素用列表，比如这里叫 Group 代表我们分析的变量因素名这一列是 Group
- `ss_type` 是平方和计算类型，一般默认2，不是我们考虑的重点，仅对不平衡设计的多因素ANOVA有效
- `detaild` 是一个bool类型变量，默认为True，输出详细表格。
- `effsize` 是效应量，默认为 np2：$\eta^2_p$，或者可以是 n2：$\eta^2$，不过对于单因素ANOVA都等于 $\frac{SS_{Between}}{SS_{Total}}$

输出的结果是表格的形式，包含以下参数：
- `Source`：变异来源，本例中为Group和Within
- `SS`：平方和
- `DF`：自由度
- `MS`：均方（SS/DF）
- `F`：F统计量，$ F=F(df1, df2) $
- `p_unc`：未校正 $p$ 值
- `np2`：$偏 \eta^2$，参考Cohen（1988）的基准：
  | 效应大小 | $\eta_p^2$ 区间          |
  |----------|--------------------------|
  | 弱效应   | $\eta_p^2 < 0.01$        |
  | 小效应   | $0.01 \le \eta_p^2 < 0.06$ |
  | 中等效应 | $0.06 \le \eta_p^2 < 0.14$ |
  | 大效应   | $\eta_p^2 \ge 0.14$      |

此时我们得到未检验的 $p$ 值，与$H_0$中的$ \alpha$ 值相比。如果$ p>\alpha$ 则不拒绝$H_0$，认为是没有差异的；若$ p<\alpha$ 则拒绝$H_0$，接下来进行事后检验两两之间差异的显著性。

### 2.Post-Hoc多重检验
在多重检验中，会出现假阳性概率叠加的情况，为了控制假阳性，我们会对两两比较的$ p $值进行校正。
最常用的是`pairwise_tukey()`，在方差齐性时使用：

```py
data.pairwise_tukey(
    dv='x',  
    between='Group',  
    alpha=0.05,      
    tail='two-sided'  
)
```
- `between`是分组变量名称
- `alpha`是显著性阈值，一般默认为0.05
- `tail`是检验类型，默认为`two-sided` 双侧检验，`greater` 和 `less `是单侧检验

输出的表格也有以下参数：
- `A`：第一组名称
- `B`：第二组名称
- `mean(A)`：A组均值
- `mean(B)`：B组均值
- `diff`：A组均值减B组均值
- `se`：标准误，$SE=\sqrt{\frac{2 \cdot MSW}{n}}$
- `T`：两两比较的t值，$T=\frac{diff}{SE} $
- `p-tukey`：tukey矫正后的p值
- `hedges`：Hedges g效应量是校正了小样本偏差的标准化均值差

## RM one-way ANOVA
### 1.数学原理
与非重复的单因素方差分析不同，对于同一个对象，我们有不同情况下的多个观测值。
对于非重复的单因素方差分析，$SST=SSB+SSE$，如果进一步讲被试均值的$SSS$（$SS_{Subject}$）分解出来，则低于重复单因素方差分析有：$SST=SSB+SSS+SSE_{rm}$

假设有n个受试者，重复测量k次，$x_{ij}$代表第i个受试者第j次实验的观测值，总观测数为N=kn：
$$SST=\sum_{i=1}^{n}\sum_{j=1}^{k}(x_{ij}-\bar{x}_{\cdot \cdot})^2, df_t=N-1$$

$$SSB=\sum_{j=1}^{k}n\cdot (\bar{x}_{\cdot j} - \bar{x}_{\cdot \cdot})^2, df_b=k-1$$

$$SSS=\sum_{i=1}^{n}k \cdot (\bar{x}_{i \cdot}-\bar{x}_{\cdot \cdot})^2, df_s=n-1$$

$$SSE_{rm}= \sum_{i=1}^{n}\sum_{j=1}^{k}[(x_{ij}-\bar{x}_{i\cdot})-(\bar{x}_{\cdot j} - \bar{x}_{\cdot \cdot})]^2, df_e=df_t-df_b-df_s=(n-1)\cdot(k-1)$$

然后对上面的平方和进行F检验：

- 组内效应：$$F(df_b, df_e)=\frac{SSB/df_b}{SSE_{rm}/df_e} $$ 

- 个体效应：$$F(df_s, df_e)=\frac{SSS/df_s}{SSE_{rm}/df_e} $$ 

一般情况下，我们都更关注组内效应。

### 2.Python函数

```py
pg.rm_anova(
  data=data, 
  dv='MoodGain', 
  within='Time', 
  subject='ID',
  detailed=True, 
  effsize='np2',
  correction='auto'
)
```
这里只需要关注新增的三个个参数：
- `within`：重复测量的条件、时间点
- `subject`：被试者编号，每个被试者的编号都唯一
- `correction`：球形性矫正，默认为 auto，自动检验球形性，不满足则用GG矫正；此外还有GG和HF两种矫正方式

而输出同样是列表：
- `Source`：变异来源（Within 因素、Error）
- `SS`：平方和
- `DF`：自由度
- `MS`：均方
- `F`：F 统计量
- `p_unc`：未校正的 p 值
- `np2`：$偏\eta^2$ 效应量
- `eps`：衡量球形度受损程度。若 $p < 0.05$ 且 $eps > 0.75$，说明球形度偏离较轻，可用 HF 校正；若 eps 较低则偏离严重，需用 GG 校正。

### 3.Post-Hoc事后检验

同样的，一样可以进行事后多重检验：
```py
pg.pairwise_tests(
  data=df_long,
  dv='Score',
  within='Training', 
  subject='Subject',
  parametric=True,
  padjust='bonf'
  )
```
这里主要介绍新增的参数：
- `parametric=True`，代表使用配对t检验
- `padjust`代表多重比较的较真方法，有 `bonf` Bonferroni矫正、`none` 不矫正、`holm` Holm矫正、`sid` Sidak矫正、`fdr_bh `FDR/BH矫正

输出列表：
- `Contrast`：对比的因素名称
- `A`：第一组/条件名称
- `B`：第二组/条件名称
- `Paired`：是否为配对/重复测量（`True` 表示重复测量，对应 RM ANOVA）
- `Parametric`：是否为参数检验（`True `表示使用 t 检验）
- `T`：t 统计量
- `dof`：自由度
- `alternative`：检验方向（two-sided 为双侧检验）
- `p-unc`：未校正的 p 值
- `p-corr`：多重比较校正后的 p 值
- `p-adjust`：多重比较校正方法
- `BF10`：Bayes 因子（可选输出项，非必需指标）
- `hedges`：Hedges' g 效应量

## 方差分析的具体步骤
### 1.ANOVA的适用条件

#### 1)数据类型与独立性
适用于连续数值数据，且组内观测值相互独立。

#### 2)数据正态性要求
一般采用以下函数来检验正态性：

  ```py
  pg.normality(
    data=data, 
    dv="Score", 
    group="Test",
    alpha=0.05
    )
  ```

#### 3)方差齐性要求
一般采用以下函数来检验正态性：
  ```py
  pg.homoscedasticity(
    data=data, 
    dv="Scores", 
    group="Group",
    method="levene" 
    )
  ```
  其中方法的选取有 `levene` 和 `bartlett` 检验两种，前者是默认首选，数据不满足正态分布也能使用。

#### 4)重复测量设计的球形性要求
  对于重复测量ANOVA，针对同⼀组观测对象不同条件下的观测值的差值⽅差相等。采用一下函数：
  ```py
  pg.sphericity(
    data=aphasia,
    dv="Score",
    within="Test",
    subject="Subject"
    ) 
  ```
  输出一下数据：
  - `spher`：是否满足球形性，True 满足，False 不满足
  - `W`：Mauchly's W 统计量，越接近 1球形性越好
  - `chi2`：卡方统计量
  - `dof`：自由度
  - `pval`：球形性检验的 p 值，一般认为大于0.05就满足球形性

#### 综上，不同条件下选取的函数推荐如下：

| 情况 | 推荐方法 | Python 函数 |
| :--- | :--- | :--- |
| 数据近似正态、方差齐性 | 经典 ANOVA | `pg.anova()` |
| 方差不齐或样本不平衡 | Welch ANOVA | `pg.welch_anova()` |
| 非正态、独立样本 | Kruskal-Wallis | `pg.kruskal()` |
| 非正态、重复测量 | Friedman 检验 | `pg.friedman()` |

### 2.结果汇报
可以参考以下APA格式进行汇报：
- 描述性统计结果
- 方差分析 F 检验结果：自由度、F统计量的值、p值、效应量等等
- 事后检验结果
- 解释统计结果的生物医学意义

## 2-way ANOVA
### 数学模型
对于没有重复测量的 2-way ANOVA：
我们设两个因素分别为A和B，各有a、b个水平，每个水平上进行n次独立观测，则一共有 $N=a\times b \times n$ 个观测值。

$$
SST=SSA+SSB+SSAB+SSE
$$

- 总平方和：
\[ SST = \sum_{i=1}^{a} \sum_{j=1}^{b} \sum_{k=1}^{r} (y_{ijk} - \bar{y}_{...})^2, df_t=N-1 \]

- 因素A平方和：
\[ SSA = r \cdot b \cdot \sum_{i=1}^{a} (\bar{y}_{i..} - \bar{y}_{...})^2 , df_a=a-1 \]

- 因素B平方和：
\[ SSB = r \cdot a \cdot \sum_{j=1}^{b} (\bar{y}_{.j.} - \bar{y}_{...})^2, df_b=b-1 \]

- 交互作用平方和：
\[ SSAB = r \cdot \sum_{i=1}^{a} \sum_{j=1}^{b} (\bar{y}_{ij.} - \bar{y}_{i..} - \bar{y}_{.j.} + \bar{y}_{...})^2 , df_{ab}=(a-1)(b-1) \]

- 误差平方和：
\[ SSE = \sum_{i=1}^{a} \sum_{j=1}^{b} \sum_{k=1}^{r} (y_{ijk} - \bar{y}_{ij.})^2, df_e=ab(r-1) \]

然后对三个因素进行F检验：
- 因素A； $$ F_a(df_a, df_e)=\frac{SSA/df_a}{SSE/df_e} $$

- 因素B： $$ F_a(df_b, df_e)=\frac{SSB/df_b}{SSE/df_e} $$

- 因素A与因素B交互作用： $$ F_{ab}(df_{ab}, df_e)=\frac{SSAB/df_{ab}}{SSE/df_e} $$

效应量：
$$
  \eta^2_p=\frac{SSA}{SSA+SSE}
$$

### 其他
2-way ANOVA 与 1-way ANOVA 的核心逻辑高度一致，但其分析的重点在于交互作用。在生物医学研究中，建议遵循以下分析路径： 
1.  交互作用优先：在查看结果时，应首先关注交互项的 $p$ 值。若交互项显著 ($p < 0.05$)，说明因素 A 的效应依赖于因素 B 的水平，此时单独解释主效应可能会造成误导，应进一步进行主效应分析，即在固定因素 B 某个水平的情况下，讨论因素 A 的差异。  

2. RM 逻辑的延展：对于 2-way RM ANOVA，其核心依然是通过数学分解剥离出“受试者个体差异” ($SSS$)。分析思路与 1-way RM 保持一致，只是变异来源中增加了第二个因素及其产生的交互项。 
