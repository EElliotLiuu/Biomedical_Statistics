
# 方差分析（ANOVA）

## 1. ANOVA 概述
方差分析（Analysis of Variance，**ANOVA**）是用于**比较两个及以上独立组均值差异**的统计方法，是t检验的扩展。
- 核心逻辑：通过分解**总变异**为**组间变异**和**组内变异**，检验变异差异是否显著
- 适用场景：因变量为连续变量，自变量为分类变量（2个及以上水平）
- 分类：单因素ANOVA、双因素ANOVA、重复测量ANOVA

---

## 2. 单因素方差分析（1-way ANOVA）
### 2.1 核心假设
1. **独立性**：各组样本相互独立
2. **正态性**：各组数据服从正态分布
3. **方差齐性**：各组总体方差相等

### 2.2 统计假设
- 原假设 $H_0$：$\mu_1 = \mu_2 = \dots = \mu_k$（所有组均值相等）
- 备择假设 $H_1$：**至少两组均值不相等**

### 2.3 核心公式（标准推导）
总平方和分解：
$$SST = SSB + SSE$$

| 指标         | 全称           | 公式                                          | 自由度  |
|--------------|----------------|-----------------------------------------------|---------|
| $SST$        | 总平方和       | $\sum_{i=1}^{k}\sum_{j=1}^{n_i}(x_{ij}-\bar{x}_{..})^2$ | $N-1$   |
| $SSB$        | 组间平方和     | $\sum_{i=1}^{k}n_i(\bar{x}_{i.}-\bar{x}_{..})^2$         | $k-1$   |
| $SSE$        | 组内/误差平方和 | $\sum_{i=1}^{k}\sum_{j=1}^{n_i}(x_{ij}-\bar{x}_{i.})^2$ | $N-k$   |

均方与F统计量：
$$MSB = \frac{SSB}{k-1}, \quad MSE = \frac{SSE}{N-k}$$
$$F = \frac{MSB}{MSE}$$

### 2.4 效应量
- 单因素ANOVA：$\eta^2 = \frac{SSB}{SST}$（偏η²与η²等价）
- 解读标准：0.01（小）、0.06（中）、0.14（大）

---

## 3. 双因素方差分析（2-way ANOVA）
### 3.1 核心概念
同时分析**两个自变量**对因变量的影响，包含三类效应：
1. 因子A的**主效应**
2. 因子B的**主效应**
3. A×B的**交互作用**（核心重点）

### 3.2 平方和分解
$$SST = SS_A + SS_B + SS_{A×B} + SSE$$

### 3.3 适用场景
研究两个分类变量的独立效应与交互效应（例：药物剂量+给药时间对疗效的影响）

---

## 4. 重复测量方差分析（RM ANOVA）
### 4.1 定义
重复测量ANOVA（Repeated Measures ANOVA）：**同一被试在不同条件/时间下的多次测量**，解决配对样本的多组比较问题。

### 4.2 核心检验：球形性检验
- 原假设：数据满足球形性假设
- 违反处理：使用**球形性矫正**
  - Greenhouse-Geisser（GG）矫正：$\varepsilon<0.75$ 时使用
  - Huynh-Feldt（HF）矫正：$0.75 \leq \varepsilon<0.9$ 时使用

### 4.3 平方和分解
$$SST = SS_{被试间} + SS_{被试内}$$
$$SS_{被试内} = SS_{处理} + SS_{误差}$$

---

## 5. Python 实战代码
### 5.1 依赖安装
```python
# 安装核心库
pip install pingouin scipy pandas numpy
```

### 5.2 基础数据格式
```python
import pandas as pd
import pingouin as pg

# 示例数据集
data = pd.DataFrame({
    "ID": range(1, 31),          # 被试编号
    "Group": ["A"]*10 + ["B"]*10 + ["C"]*10,  # 分组（3组）
    "Score": [2.1,2.3,1.9,2.5,2.2,3.0,2.8,2.4,2.6,2.0,
              3.2,3.5,3.1,3.3,3.0,3.8,3.6,3.4,3.2,3.7,
              4.1,4.3,4.0,4.2,4.5,4.8,4.6,4.4,4.2,4.7]
})
```

### 5.3 单因素ANOVA 完整流程
```python
# 1. 方差齐性检验
levene = pg.homoscedasticity(data, dv="Score", group="Group")
print("方差齐性检验：\n", levene)

# 2. 单因素ANOVA
anova_res = pg.anova(data=data, dv="Score", between="Group", detailed=True, effsize="np2")
print("\n单因素ANOVA结果：\n", anova_res)

# 3. 事后检验（Tukey法，多重比较校正）
post_hoc = pg.pairwise_tukey(data=data, dv="Score", between="Group")
print("\n事后检验结果：\n", post_hoc)
```

### 5.4 重复测量ANOVA
```python
# 构造重复测量数据（同一被试3次测量）
rm_data = pd.DataFrame({
    "ID": list(range(1,11))*3,
    "Time": ["T1"]*10 + ["T2"]*10 + ["T3"]*10,
    "Score": [1.5,1.6,1.4,1.7,1.5,1.8,1.6,1.5,1.7,1.6,
              2.2,2.3,2.1,2.4,2.2,2.5,2.3,2.2,2.4,2.3,
              2.8,2.9,2.7,3.0,2.8,3.1,2.9,2.8,3.0,2.9]
})

# 重复测量ANOVA + 球形性检验 + 矫正
rm_anova = pg.rm_anova(data=rm_data, dv="Score", within="Time", subject="ID", detailed=True)
print("重复测量ANOVA结果：\n", rm_anova)
```

---

## 6. 结果解读与规范汇报
### 6.1 ANOVA 结果核心指标
- $F$：F统计量
- $p$：p值（$p<0.05$ 表示差异显著）
- $\eta_p^2$：偏效应量（效应大小）
- 自由度：$(df_{组间}, df_{误差})$

### 6.2 APA 规范汇报示例
> 单因素ANOVA结果显示，三组得分差异显著，$F(2, 27) = 158.2, p < 0.001, \eta_p^2 = 0.92$。事后Tukey检验表明，A组、B组、C组两两之间差异均显著（$p<0.001$）。

---

## 7. 常见问题与解决方案
1. **方差齐性违反**
   - 替代方法：Welch ANOVA
2. **正态性违反**
   - 替代方法：Kruskal-Wallis H检验（非参数）
3. **ANOVA显著但事后检验不显著**
   - 原因：样本量小、多重比较校正保守
   - 方案：更换校正方法（Holm/FDR-BH）
4. **重复测量球形性违反**
   - 方案：使用GG/HF矫正结果

---

## 附录：术语对照表
| 英文缩写 | 全称                     | 中文释义       |
|----------|--------------------------|----------------|
| ANOVA    | Analysis of Variance    | 方差分析       |
| SSB      | Sum of Squares Between   | 组间平方和     |
| SSE      | Sum of Squares Error    | 误差平方和     |
| RM       | Repeated Measures        | 重复测量       |
| GG       | Greenhouse-Geisser       | 球形性矫正方法 |



