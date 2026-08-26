# DropoutKiller：已知 dropout 坐标后的局部表达恢复数学架构

> 状态：2026-08-26。本文只讨论 **dropout 坐标已经确定以后** 的表达恢复问题。检测不是本文件的优化目标。

## 1. 问题重新定义

设归一化后的基因 × 细胞表达矩阵为

\[
X\in\mathbb R_+^{G\times N},
\]

已知技术性 dropout 坐标集合

\[
\mathcal M=\{(g,c):X_{gc}=0,\;(g,c)\text{ 已被判定为 technical dropout}\}.
\]

同时已经有：

- 最终 SuperCell membership：\(m_c\)；
- biological embedding：\(z_c\in\mathbb R^d\)；
- SuperCell / walktrap hierarchy：\(\mathcal T\)；
- 可选 hard biological stratum：\(s_c\)，例如 major cell type、condition、donor 等；
- 非 dropout 坐标均保持原值，不允许恢复器覆盖。

因此真正需要估计的是

\[
\widehat X_{gc},\qquad (g,c)\in\mathcal M,
\]

而不是重新学习整个表达矩阵，也不需要恢复器再次判别 biological zero 与 dropout。

在当前问题中，一个自然 estimand 是

\[
E[X_{gc}\mid (g,c)\in\mathcal M,\;z_c,\;m_c,\;\mathcal T,\;s_c,\;\text{reliable observations}],
\]

若 dropout 判定已经意味着“该坐标原本应有可检测表达”，则更保守的 positive-conditional estimand 为

\[
E[X_{gc}\mid X_{gc}>0,\;z_c,\;m_c,\;\mathcal T,\;s_c,\;\text{reliable observations}].
\]

这与大部分经典方法面对的原始问题不同：经典方法通常同时解决 dropout detection、global manifold learning 和 imputation；这里前两部分已经由外部结构提供。

---

## 2. 已发表和近期预印本方法的核心数学母架构

下面不是简单按软件罗列，而是按可迁移的数学结构归类。

### 2.1 Cluster / neighborhood smoothing

#### DrImpute

先重复聚类，再用目标细胞所在簇内同一基因的平均表达恢复：

\[
\widehat X_{gc}=\frac{1}{|C(c)|}\sum_{j\in C(c)}X_{gj}.
\]

优点是稳定、快速；缺点是 membership 内所有细胞近似 exchangeable，容易压低 differential variability。

**对当前问题的价值：** 适合作为最简单 oracle baseline，而不适合作为最终默认模型。

#### MAGIC

构建细胞 affinity，再归一化得到 Markov 转移矩阵 \(P\)，做 \(t\) 步扩散：

\[
\widehat X=P^tX.
\]

其本质是 graph heat diffusion / low-pass filtering。

**对当前问题的价值：** 可以借用“图上的局部传播”思想，但不应直接对整个矩阵扩散，因为会修改可靠非 dropout 坐标并造成过平滑。

参考：van Dijk et al., *Cell*, 2018, doi:10.1016/j.cell.2018.05.061.

### 2.2 Cell self-representation：最接近当前问题的一类

#### scImpute

scImpute 先用 Gamma-Normal mixture 得到 dropout probability；在恢复阶段，对于目标细胞 \(c\)，只用可信基因集合 \(B_c\) 学其他细胞的非负表示权重：

\[
\widehat\beta^{(c)}=
\arg\min_{\beta\ge0}
\left\|X_{B_c,c}-X_{B_c,N_c}\beta\right\|_2^2.
\]

然后对需要插补的基因使用同一组 cell weights：

\[
\widehat X_{gc}=X_{g,N_c}\widehat\beta^{(c)}.
\]

这是一个非常重要的设计：**目标基因本身不需要参与 donor weight 学习；cell-specific weights 可以复用于该细胞的多个 dropout 基因。**

参考：Li & Li, *Nature Communications*, 2018, doi:10.1038/s41467-018-03405-7.

#### VIPER

VIPER 将一个细胞表示为少量邻居细胞的非负凸组合：

\[
X_{\cdot c}\approx \sum_{j\ne c}b_{cj}X_{\cdot j},
\]

约束

\[
b_{cj}\ge0,\qquad \sum_j b_{cj}=1,
\]

并通过候选预筛选 + 稀疏非负回归得到少量 donor。对缺失基因：

\[
\widehat X_{gc}=\sum_j b_{cj}X_{gj}.
\]

VIPER 的核心优点是 cell-specific、asymmetric、sparse donor representation，并明确以 preserving variability 为设计目标。

参考：Chen & Zhou, *Genome Biology*, 2018, doi:10.1186/s13059-018-1575-1.

#### SH-CoCL 的最终局部恢复器，2026

SH-CoCL 前半段用 stochastic hypergraph + contrastive learning 学 embedding，但其最终 imputation 阶段本身非常值得单独拆出来。对 cell \(c\)，先在 embedding 中寻找 \(k\) 个邻居。令 \(X_c^{obs}\) 是这些邻居在 query cell 可靠 observed genes 上的表达矩阵，\(y_c\) 是 query cell 对应 observed expression，则权重由 ridge self-reconstruction 得到：

\[
\widehat w_c=
\left[(X_c^{obs})^TX_c^{obs}+\lambda I\right]^{-1}(X_c^{obs})^Ty_c.
\]

然后把**同一组权重**直接作用到 missing genes 的邻居表达矩阵 \(X_c^{mis}\)：

\[
\widehat x_c^{mis}=X_c^{mis}\widehat w_c.
\]

这与 scImpute/VIPER 属于同一个非常关键的数学母架构：**在可靠坐标上学习 cell-specific donor weights，再把权重迁移到已知 missing/dropout 坐标。**

SH-CoCL 为了得到 embedding 还训练概率超图：

\[
p_{ij}=\sigma(s_{ij}),\qquad
\widetilde h_{ij}=\sigma\left(\frac{s_{ij}+\epsilon_{ij}}{\tau}\right),
\]

并通过 hypergraph convolution / contrastive learning 得到 \(z_c\)。但在当前 DropoutKiller 设定中，biological embedding 与 SuperCell hierarchy 已经存在，因此无需重训这一大段；可以直接保留其最后的 local self-representation recovery。

参考：Zi et al., *Complex & Intelligent Systems*, 2026, doi:10.1007/s40747-026-02360-x.

### 2.3 Bayesian shrinkage / posterior recovery

#### SAVER

对 UMI count：

\[
Y_{gc}\mid\lambda_{gc}\sim\operatorname{Poisson}(s_c\lambda_{gc}),
\]

\[
\lambda_{gc}\sim\operatorname{Gamma}(\alpha_{gc},\beta_{gc}).
\]

先用其他基因的 Poisson-Lasso 预测 prior mean \(\mu_{gc}\)，再得到 Gamma posterior。后验均值可写成 observed normalized count 与 prediction 的精度加权：

\[
E[\lambda_{gc}\mid Y]
=w_{gc}\frac{Y_{gc}}{s_c}+(1-w_{gc})\mu_{gc}.
\]

SAVER 的重要贡献不是“均值插补”，而是**保留 posterior uncertainty**。

参考：Huang et al., *Nature Methods*, 2018, doi:10.1038/s41592-018-0033-z.

#### SAVER-X

用 transfer-learning autoencoder 给出 \(\widehat\Lambda_{gc}\)，再做 empirical-Bayes shrinkage：

\[
\widehat X_{gc}=
\frac{\sigma_{2g}^2}{\sigma_{1g}^2+\sigma_{2g}^2}\frac{Y_{gc}}{l_c}
+
\frac{\sigma_{1g}^2}{\sigma_{1g}^2+\sigma_{2g}^2}\widehat\Lambda_{gc}.
\]

**对当前问题的价值：** 不必照搬 count model，但必须保留“预测均值 + 预测方差/后验抽样”的思想，否则 deterministic imputation 天然压缩 DV。

### 2.4 Low-rank / matrix-factorization / sparse recovery

#### ALRA

先做 rank-\(k\) truncated SVD：

\[
X\approx U_k\Sigma_kV_k^T=\widetilde X,
\]

再利用低秩近似中负值的分布，对每个基因自适应阈值，把接近零的重构值重新压回零。

其关键假设是：真实矩阵近似 low-rank、非负，并含 biological zeros。

**对当前问题的价值：** 若 membership 已经把细胞限制到局部状态空间，可以把 low-rank 从“全局假设”变成“局部候选 engine”，风险显著降低；当前 `masked_factor` / `tree_local_factor` 已覆盖这一方向的一部分。

参考：Linderman et al., *Nature Communications*, 2022, doi:10.1038/s41467-021-27729-z.

#### scIRT，2026

scIRT 先融合多种距离的 rank：

\[
M=\sum_n t_nR^{(n)},
\]

再用与相似细胞的 SMOTE 型随机插值：

\[
a_i^*=a_i+\sum_j f_{ij}(a_j-a_i),
\]

其中 donor 权重由相似性 rank 和几何分布随机化；随后用 NMF 反复约束：

\[
A^+\approx Z^+H^+,
\qquad
\min_{Z,H\ge0}\|A-ZH\|_F^2.
\]

其价值在于“多尺度 neighborhood borrowing + low-rank regularization”可以分层执行；缺点是随机插值和 NMF 都可能进一步平滑，因此不宜直接作为 DV-sensitive 默认恢复器。

参考：Mou et al., *International Journal of Molecular Sciences*, 2026, doi:10.3390/ijms27031173.

#### GSCI，2026

GSCI 把生成式预恢复与 sparse compressed sensing 串联。其稀疏层可抽象为：

\[
x\approx D w,
\]

\[
\widehat w=\arg\min_w\|x-Dw\|_2^2+\lambda\|w\|_1,
\]

再由固定/学习 dictionary 重构 \(\widehat x=D\widehat w\)。完整方法先用共享 encoder + VAE/GAN 得到初始恢复，再通过 SVD dictionary、LASSO/OMP 和 compressed-sensing optimization 做 sparse refinement，并使用 zero-aware reconstruction objective。

**对当前问题的价值：** sparse representation 值得保留，但全局 VAE/GAN 训练并非必要；更自然的局部版本是只在 biological neighborhood 内建立 small dictionary，并仅对 oracle dropout coordinates 读取恢复值。

参考：Tian et al., *Computational and Structural Biotechnology Journal*, 2026, doi:10.34133/csbj.0124.

### 2.5 Two-sided self-representation

#### scTSSR

cell-side representation：

\[
X_{gc}=\sum_{c'\ne c}X_{gc'}B_{c'c}+E_{gc}.
\]

gene-side representation：

\[
X_{gc}=\sum_{g'\ne g}A_{gg'}X_{g'c}+E_{gc}.
\]

两者组合：

\[
X\approx AX+XB+AXB.
\]

并用稀疏约束学习 \(A,B\)。

**对当前问题的价值：** 说明“cell similarity”和“gene coexpression”可以是两个独立专家，而不是必须融合成一个 latent factor。当前问题中已有可靠 biological embedding，因此 cell-side 部分尤其值得保留；gene-side 可以作为 residual expert，而不应强迫进入主恢复器。

参考：Gong et al., *Bioinformatics*, 2020, doi:10.1093/bioinformatics/btaa108.

### 2.6 Multi-structure / mixture-of-experts

#### scMOO

同时学习 gene-side structure \(A\)、cell-side structure \(B\) 和 low-rank structure \(L\)：

\[
f_1=\|X-AY\|_F^2+\lambda_1\|A\|_1,
\]

\[
f_2=\|X-YB\|_F^2+\lambda_1\|B\|_1,
\]

\[
f_3=\|X-L\|_F^2.
\]

总体目标通过 \(\omega_1,\omega_2,\omega_3\ge0\)、\(\sum\omega_i=1\) 自动组合多个结构。其 element-wise 更新可以解释为 observed、gene-side、cell-side、low-rank 四者的加权平均。

**对当前问题的价值：** 非常适合转化为“多个局部恢复器的自监督 stacking”，而不是在全矩阵上联合求解所有结构。

参考：Gong et al., *Bioinformatics*, 2022, doi:10.1093/bioinformatics/btac300.

### 2.7 Deep count / graph / autoencoder family

#### DCA

decoder 输出 ZINB 参数：

\[
X_{gc}\sim\operatorname{ZINB}(\pi_{gc},\mu_{gc},\theta_{gc}),
\]

训练目标为 ZINB negative log-likelihood；\(\mu_{gc}\) 作为 denoised expression。

#### DeepImpute

把 target genes 分块，每个子网络只用与 target genes 高相关、且不属于 target set 的 genes 作为输入：

\[
\widehat X_{T,c}=f_\theta(X_{P,c}),\qquad P\cap T=\varnothing.
\]

这与当前“目标泄漏必须避免”的要求高度一致。

#### scGNN / GraphSCI

典型 graph convolution：

\[
H^{(l+1)}=\sigma(\widetilde D^{-1/2}\widetilde A\widetilde D^{-1/2}H^{(l)}W^{(l)}),
\]

再由 autoencoder / graph-regularized decoder 重构表达。

#### scGImpute，2026

scGImpute 属于 KNN + 双层 multi-head graph attention 的 hybrid 路径。标准 GAT 可抽象为

\[
e_{cj}^{(h)}=a_h(W_hz_c,W_hz_j),
\qquad
\alpha_{cj}^{(h)}=\operatorname{softmax}_{j\in N(c)}e_{cj}^{(h)},
\]

\[
z_c'=\Vert_h\sigma\left(\sum_{j\in N(c)}\alpha_{cj}^{(h)}W_hz_j\right).
\]

其核心价值不是必须引入神经网络，而是说明 fixed-distance KNN weights 可以进一步变成 **data-adaptive attention weights**。在当前问题里，可把 attention 限制在既有 biological neighborhood 内并做轻量化，而无需重新学习全局图。

参考：Borah & Das, *Computational Biology and Chemistry*, 2026, doi:10.1016/j.compbiolchem.2025.108856.

**对当前问题的总体判断：** 如果 biological embedding 与 SuperCell hierarchy 已经可靠，重新训练 global GNN/VAE 往往是在重复学习已经存在的信息；深模型更适合作为 comparator，或者只保留其“attention / uncertainty / nonlinear local predictor”模块。

### 2.8 2025–2026 selective / local recovery 方法

#### scTsI，2025，已发表

第一阶段用邻近 cell 和邻近 gene 均值：

\[
X_{ij}^{(1)}=
\frac12\left(
\frac1{k_1}\sum_{u=1}^{k_1}X_{iu}
+
\frac1{k_2}\sum_{u=1}^{k_2}X_{uj}
\right).
\]

第二阶段用 bulk RNA 约束：

\[
X^{(2)}=\arg\min_X
\|X-X^{(1)}\|_F^2+\lambda\|Xa-d\|_2^2.
\]

它明确只更新 zero entries。

对当前问题而言，第一阶段的 cell-local borrowing 可保留；bulk constraint 不是必需条件。

参考：Zhang et al., *Briefings in Bioinformatics*, 2025, doi:10.1093/bib/bbaf298.

#### D3Impute，2025，已发表

三部分：distribution-aware normalization、dropout discriminator、density-guided local imputation。当前问题已经完成 dropout discrimination，因此最有价值的是其思想：**恢复阶段应该保持局部 density / neighborhood，而不是做全局平滑。**

参考：Huang et al., *PLOS Computational Biology*, 2025, doi:10.1371/journal.pcbi.1013744.

#### scZN，2026，已发表

以 transcription burst + dropout 为 ZINB 统计层，同时做 biologically constrained NMF：

\[
X\approx WH,\qquad W,H\ge0,
\]

并叠加 ZINB NLL、NMF reconstruction、variance alignment、classification/prior regularization。

对当前问题而言，global NMF/ZINB 过重，但“已有 biological prior 可以作为 soft constraint，而不是 hard label”这一点值得吸收。

参考：Wu et al., *PLOS Computational Biology*, 2026, doi:10.1371/journal.pcbi.1014051.

#### scZiva，2026，已发表

VAE decoder：

\[
X_{gc}\sim\operatorname{ZINB}(\pi_{gc},\mu_{gc},\theta_{gc}),
\]

总损失：

\[
\mathcal L=\mathcal L_{ZINB}+\mathcal L_{KL}+\lambda\mathcal L_{MSE},
\]

其中 MSE 只施加在 observed coordinates。其 selective recovery 为

\[
\widehat X_{gc}=(1-\pi_{gc})\mu_{gc}
\]

仅在零值且 dropout probability 超阈值时使用。

参考：Vo et al., *BMC Bioinformatics*, 2026, doi:10.1186/s12859-026-06422-2.

#### scDDI，2026，已发表

先 PNB dropout detection，再计算 dropout-aware weighted cosine similarity，选择局部细胞，最后用 decision-tree regression 做局部非线性恢复。

当前问题已经不需要 PNB detection，因此可以直接吸收：

\[
\widehat X_{gc}=f_{g,c}^{tree}(X_{D(c),\,-g}),
\]

其中训练样本只来自局部相似细胞。

参考：*Briefings in Bioinformatics*, 2026, bbag072.

#### SH-CoCL，2026，已发表

SH-CoCL 用 Gumbel-Sigmoid 概率超图和 co-contrastive learning 学稳健 embedding：

\[
\widetilde h_{ij}=\sigma\left(\frac{s_{ij}+\epsilon_{ij}}{\tau}\right),
\]

超图卷积为

\[
X^{(l+1)}=\sigma\left(D_v^{-1/2}\widetilde H D_e^{-1}\widetilde H^TD_v^{-1/2}X^{(l)}W^{(l)}\right).
\]

但对当前问题最关键的是：其最终表达恢复仍然回到**邻居 + observed-gene ridge self-representation + 同权重 missing-gene prediction**，因此为当前局部 recovery-only 设计提供了很直接的近期证据。

参考：Zi et al., *Complex & Intelligent Systems*, 2026, doi:10.1007/s40747-026-02360-x.

#### GSCI，2026，已发表

GSCI 是 generative pre-recovery + sparse compressed-sensing refinement 的两阶段框架。它提示另一个可探索方向：在已经有局部 biological neighborhood 的情况下，可省略 global generator，只在邻域中建立稀疏 dictionary，解

\[
\widehat w=\arg\min_w\|P_R(x_c-D_cw)\|_2^2+\lambda\|w\|_1,
\]

然后只对 dropout coordinates 读取 \(D_c\widehat w\)。

参考：Tian et al., *Computational and Structural Biotechnology Journal*, 2026, doi:10.34133/csbj.0124.

#### SCR-MF，2025，预印本

把 scRecover 的 dropout detection 与 missForest / Random Forest imputation 分开。对当前问题的意义主要是证明“**先确定缺失坐标，再用非参数局部 predictor**”是合理的模块边界。

预印本：arXiv:2511.16923.

#### 2026 large-scale benchmark，预印本

30 个数据集、10 种实验 protocol、15 个 imputation method、6 类 downstream task 的比较显示：

- traditional model/smoothing/low-rank methods总体并不弱于、且常优于 DL；
- numerical recovery 最优不等于 biological downstream 最优；
- performance 强烈依赖 dataset、protocol 和 task；
- 没有一个方法在所有场景持续最优。

因此工程设计不应该预设单一算法，而应保留可比较的多个 recovery engine 和统一 oracle-mask benchmark。

预印本：Iwashita et al., arXiv:2603.24626.

---

## 3. 当前约束下可以进一步拓展的数学架构

## A. Nadaraya–Watson local kernel mean

先定义 cell geometry kernel：

\[
k_{cj}=\exp\left[-
\alpha_E\frac{d_E(c,j)^2}{2h_c^2}
-\alpha_T\frac{d_T(c,j)}{\tau_T}
-\eta I(m_c\ne m_j)
\right].
\]

对 dropout target 只保留 positive donors：

\[
D_{gc}=\{j:X_{gj}>0,\;s_j=s_c\}.
\]

预测：

\[
\widehat X_{gc}=
\frac{\sum_{j\in D_{gc}}k_{cj}X_{gj}}
{\sum_{j\in D_{gc}}k_{cj}}.
\]

优点：简单、稳定、易 batch；缺点：局部常数模型仍会 bias toward mean。

## B. Soft-hierarchy barycentric / LLE recovery

这是 0.7.0 新增的主要 candidate engine。

### B1. hard boundary 与 soft membership

hard stratum 仍然是绝对边界：

\[
s_j\ne s_c\Rightarrow w_{cj}=0.
\]

但 membership 不再是墙，而是 penalty：

\[
I(m_j\ne m_c)\text{ 只增加 donor cost，不直接置零。}
\]

### B2. geometry prior

候选细胞先由 biological embedding 的 local kNN 给出。定义

\[
C_{cj}=
\begin{cases}
\alpha_T d_T(c,j)/\tau_T
+(1-\alpha_T)d_E(c,j)^2/(2h_c^2)
+\eta I(m_j\ne m_c), & d_T\text{ available},\\
 d_E(c,j)^2/(2h_c^2)+\eta I(m_j\ne m_c), & \text{otherwise}.
\end{cases}
\]

几何先验为

\[
p_{cj}=\frac{e^{-C_{cj}}}{\sum_l e^{-C_{cl}}}.
\]

### B3. barycentric state reconstruction

令

\[
\Delta_{cj}=\frac{z_j-z_c}{h_c}.
\]

学习 query-specific donor weights：

\[
\widehat w_c=
\arg\min_{w\ge0,\;\mathbf1^Tw=1}
\left\|\sum_jw_{cj}\Delta_{cj}\right\|_2^2
+\lambda_B\|w_c-p_c\|_2^2.
\]

第一项要求 donor convex combination 在 embedding 中重构 query cell；第二项阻止权重偏离 hierarchy / membership / distance 提供的 biological prior。

该问题为 convex quadratic program；当前实现使用 projected gradient，每一步精确投影到 probability simplex。

### B4. target-gene recovery

目标基因完全不进入 \(w_c\) 的学习。对 positive donors 重新归一化：

\[
\widetilde w_{gcj}=
\frac{w_{cj}I(X_{gj}>0)}
{\sum_l w_{cl}I(X_{gl}>0)}.
\]

最终：

\[
\widehat X_{gc}=\sum_j\widetilde w_{gcj}X_{gj}.
\]

因此一个 cell 的多个 dropout genes 共用同一个 state weight vector，只在最后一步按目标基因 positive support 重归一化。

### B5. effective donors 与 uncertainty

\[
n_{eff,gc}=\frac{1}{\sum_j\widetilde w_{gcj}^2}.
\]

加权生物变异：

\[
s_{gc}^2=
\frac{\sum_j\widetilde w_{gcj}(X_{gj}-\widehat X_{gc})^2}
{1-\sum_j\widetilde w_{gcj}^2}.
\]

working predictive variance：

\[
V_{gc}=s_{gc}^2\left(1+\frac1{n_{eff,gc}}\right).
\]

它把 donor biological spread 保留下来，而不是只输出 local mean uncertainty。

## C. Reliable-gene local self-representation

这是 scImpute、VIPER 和 2026 SH-CoCL 最终恢复阶段对当前问题最直接的改写。给定 embedding/hierarchy 先筛小 donor set \(D_c\)，再只用 query cell 的 reliable coordinates \(R_c\) 学 cell weights：

\[
\widehat w_c=
\arg\min_{w\ge0,\mathbf1^Tw=1}
\|P_{R_c}(X_{\cdot c}-X_{\cdot D_c}w)\|_2^2
+\lambda_E\|w-p_c\|_2^2.
\]

或者去掉非负/和为 1 约束，使用 SH-CoCL 型 ridge closed form：

\[
\widehat w_c=(X_R^TX_R+\lambda I)^{-1}X_R^Ty_c.
\]

然后对 dropout target：

\[
\widehat X_{gc}=X_{g,D_c}\widehat w_c.
\]

关键是 \(g\notin R_c\)，并且所有已知 dropout target 坐标都从 weight-learning feature 中屏蔽。

**理论优势：** 比 embedding-only barycentric 更直接地回答“哪些邻居能重构这个 cell 的可靠转录状态”。

**工程风险：** donor-specific missing coordinates、feature scaling、大量可靠基因上的稀疏矩阵计算和潜在 library-size confounding 都必须处理。因此适合作为下一 candidate engine，与 barycentric 在相同 oracle masks 上正面对比，而不是未经验证直接替换。

## D. Local linear / local polynomial regression

在 query-local positive donors 上拟合：

\[
(\hat\beta_0,\hat\beta)=
\arg\min_{\beta_0,\beta}
\sum_{j\in D_{gc}}w_{cj}
[X_{gj}-\beta_0-\beta^T(z_j-z_c)]^2
+\lambda\|\beta\|_2^2.
\]

因为在 query 点 \(z_j-z_c=0\)，预测直接为

\[
\widehat X_{gc}=\max(0,\hat\beta_0).
\]

比 kernel mean 更能降低 manifold slope / boundary bias，但每个 gene-query 都要做局部回归，工程成本高于 barycentric mean。

## E. Graph harmonic extension

将 dropout 坐标视为 graph 上真正的 missing node labels。对每个 gene：

\[
\min_y\frac12\sum_{(c,j)}w_{cj}(y_c-y_j)^2
\]

并固定 reliable observations \(y_O=X_{gO}\)。若图 Laplacian 分块为

\[
L=\begin{bmatrix}L_{MM}&L_{MO}\\L_{OM}&L_{OO}\end{bmatrix},
\]

则 missing 部分满足

\[
\widehat y_M=-L_{MM}^{-1}L_{MO}y_O.
\]

优点：严格利用整张局部图；缺点：Laplacian harmonic solution 天然平滑，容易压缩 DV。

## F. Graph total variation / trend filtering

用 TV 代替 quadratic Laplacian：

\[
\min_y
\frac12\|P_O(y-X_g)\|_2^2
+\lambda\sum_{(c,j)}w_{cj}|y_c-y_j|.
\]

比 harmonic extension 更能保留 sharp state boundaries 和 piecewise structure，但求解器复杂度更高。

## G. Local Gaussian process

对 embedding 上的 target gene function 建局部 GP：

\[
f_g(z)\sim GP(0,k(z,z')).
\]

对 local positive donors \(D\)：

\[
\widehat X_{gc}=k_{cD}(K_{DD}+\sigma^2I)^{-1}X_{gD},
\]

\[
V_{gc}=k_{cc}-k_{cD}(K_{DD}+\sigma^2I)^{-1}k_{Dc}.
\]

优点：天然 uncertainty、局部非线性；缺点：每个 gene-local block 都要矩阵求逆，不适合百万级 event 直接默认启用。

## H. Local masked low-rank / matrix completion

在每个 biological neighborhood \(\mathcal N_c\) 上：

\[
\min_{U,V}
\|P_\Omega(X-UV^T)\|_F^2
+\lambda(\|U\|_F^2+\|V\|_F^2).
\]

只在 known dropout mask 上读出 \(UV^T\)。当前 `masked_factor` 和 `tree_local_factor` 已属于这一族的局部化版本。

## I. Local sparse dictionary / compressed sensing

在 embedding / hierarchy 给出的局部 donor matrix \(D_c\) 上，用 reliable coordinates 求 sparse coefficient：

\[
\widehat a_c=\arg\min_a
\|P_{R_c}(x_c-D_ca)\|_2^2+\lambda\|a\|_1.
\]

只在 dropout coordinates 上输出

\[
\widehat x_{c,\mathcal M}=D_{c,\mathcal M}\widehat a_c.
\]

它是 GSCI compressed-sensing 思路的 recovery-only、local 版本；与 C 的区别主要在于显式 L1 sparsity，可自动挑少量 donor / dictionary atom。

## J. Local nonlinear tree / random forest

对应 scDDI / SCR-MF：

\[
\widehat X_{gc}=f_g(X_{c,R_c};\;D_c),
\]

其中 \(f_g\) 是只在 local donor cells 中训练的 decision tree / random forest。

优点：可捕获 nonlinear interactions；缺点：membership 规模约 100–200 时容易高方差，并且每 gene 或 gene-block 训练树的成本较高。

## K. Local attention / graph attention

只在 biological candidate set 上学习 attention：

\[
e_{cj}=a_\theta(z_c,z_j,\phi_{cj}),
\qquad
\alpha_{cj}=\frac{\exp(e_{cj})}{\sum_{l\in D_c}\exp(e_{cl})},
\]

其中 \(\phi_{cj}\) 可包含 tree distance、same-membership indicator、reliable-gene cosine similarity、local density 等。恢复为

\[
\widehat X_{gc}=\frac{\sum_j\alpha_{cj}I(X_{gj}>0)X_{gj}}
{\sum_j\alpha_{cj}I(X_{gj}>0)}.
\]

优点是能非线性学习多个相似性来源；缺点是需要独立的 self-supervised oracle-mask objective，否则 attention 很容易学到不稳定的 dataset-specific shortcut。

## L. Distribution-preserving donor resampling / local quantile recovery

如果下游关心 differential variability，单一 conditional mean 不够。可直接估计 weighted empirical CDF：

\[
\widehat F_{gc}(x)=\sum_j\widetilde w_{gcj}I(X_{gj}\le x).
\]

然后：

\[
X_{gc}^{(b)}\sim\widehat F_{gc}.
\]

也可以平滑成 Gamma / log-normal / hurdle distribution。该方向不改变 deterministic mean engine，但为 DV-sensitive analysis 提供 multiple imputation。

---

## 4. 更合适的整体工程架构

建议把整个包拆成六层，而不是把 detection 和 recovery 写成一个不可替换的大模型。

```text
Expression X
   + known dropout mask M
   + biological embedding Z
   + membership m
   + SuperCell hierarchy T
   + hard strata s
          │
          ▼
1. Candidate graph builder
   hard stratum = absolute boundary
   embedding kNN = local candidate set
   membership = soft prior
   hierarchy = multiscale prior
          │
          ▼
2. Cell-state weight learner
   kernel / barycentric / reliable-gene self-expression
   sparse dictionary / local attention
          │
          ▼
3. Target-gene conditional estimator
   positive weighted mean / local linear / GP
   factor residual / tree / forest
          │
          ▼
4. Uncertainty / distribution layer
   weighted variance / Bayesian posterior
   empirical donor CDF / multiple imputation
          │
          ▼
5. Selective writer
   only M==TRUE coordinates can change
          │
          ▼
6. Oracle-mask benchmark / engine selection / stacking
```

### 统一 recovery-engine contract

每个 engine 至少返回：

```text
prediction
prediction_sd / predictive_variance
n_donors
n_eff
local prevalence
embedding distance summary
tree distance summary
engine-specific predictability
```

并必须满足：

1. observed coordinates bitwise unchanged；
2. hard strata 不跨界；
3. target gene 不得泄漏到 donor-weight training；
4. dropout mask 为空时输出等于输入；
5. 无足够 donor support 时宁可 unavailable，也不强行生成表达；
6. donor geometry 与 target-gene magnitude estimator 尽可能解耦，便于消融与审计；
7. deterministic mean 与 uncertainty / multiple-imputation layer 分离。

---

## 5. 不预设赢家：推荐 benchmark 设计

### 5.1 Oracle positive masking

从原本可靠的 positive coordinates 中抽取一部分作为人工 missing：

\[
\mathcal M_{oracle}\subset\{(g,c):X_{gc}>0\}.
\]

隐藏后恢复，再与真实值比较。这是最直接的 recovery-only benchmark。

必须分层抽样：

- original count / expression abundance；
- gene prevalence；
- membership size；
- tree distance；
- local density；
- cell state；
- condition / donor；
- high-DV vs low-DV genes。

### 5.2 Binomial thinning

如果有 raw UMI count \(C_{gc}\)：

\[
C'_{gc}\sim\operatorname{Binomial}(C_{gc},q_c).
\]

保留 full-depth count 作为近似 reference，重新执行 normalization / embedding / membership 后评估 recovery，避免只在固定 embedding 上做过于理想化的 benchmark。

### 5.3 评价指标不能只用 RMSE

数值误差：

\[
RMSE,\ MAE,\ bias,\ Pearson,\ Spearman,\ CCC.
\]

分布与 DV：

\[
\frac{\operatorname{Var}(\widehat X_g)}{\operatorname{Var}(X_g)},
\qquad
W_1(\widehat F_g,F_g),
\]

以及 quantile error、Gini error。

生物结构：

- gene-gene correlation distortion；
- marker leakage；
- DE logFC bias；
- differential-variability statistic bias；
- trajectory / local manifold distortion；
- rare-state recall。

uncertainty：

- 50/80/95% interval coverage；
- interval width；
- calibration by n_eff / expression abundance。

### 5.4 Engine selection / stacking

不要固定全数据唯一 engine。可以在 oracle masks 上学非负 stacking weight：

\[
\widehat X_{gc}
=\sum_{r=1}^R\omega_{r,h(gc)}\widehat X_{gc}^{(r)},
\]

\[
\omega_r\ge0,\qquad\sum_r\omega_r=1,
\]

其中 \(h(gc)\) 可以是 gene prevalence × local density × n_eff 的 coarse stratum。

这相当于把 scMOO 的 multi-structure 思想改造成 recovery-only mixture-of-experts，避免重新估计整个表达矩阵。

### 5.5 为什么一定要保留多个 engine

2026 大规模 benchmark 的关键信息不是“传统方法一定最好”，而是：**恢复值误差、聚类、DE、轨迹和其他生物任务对插补偏差的敏感方向不同。** 因此同一个 engine 即使 RMSE 最优，也可能不是 DV 或 DE 最优。对 DropoutKiller 更合理的是：

```text
固定同一 oracle mask
        ↓
所有 engine 同场恢复
        ↓
按 abundance / prevalence / n_eff / biological stratum 分层评分
        ↓
决定全局 default、分层 selector 或 stacking
```

而不是先通过理论偏好把某个 engine 写死。

---

## 6. 当前 0.7.0 实现选择

本次没有把 `barycentric` 直接升级为默认值，而是作为新的 benchmark engine 加入：

```r
recovery_method = "barycentric"
```

核心原因：当前文献与 2026 大规模 benchmark 都不支持“某一种数学结构在所有数据上必然最优”。

当前 engine 的设计原则是：

1. hard biological stratum 仍为绝对边界；
2. membership 从 hard block 改为 soft penalty；
3. embedding kNN 决定局部候选；
4. SuperCell hierarchy 调节跨 membership borrowing；
5. donor weights 通过 barycentric state reconstruction 学习；
6. target gene 不参与 weight learning；
7. target expression 只由 reliable positive donors 得到；
8. 同一 cell 的所有 dropout genes 共用 cell-state weights；
9. 输出 local variance、Kish n_eff 和 predictive SD；
10. 原有 `tree_local_factor`、`masked_factor`、`neighbor` 保留作为统一 benchmark comparator。

### 推荐首先比较的恢复器

1. positive membership mean；
2. embedding-only Gaussian kernel mean；
3. tree + embedding kernel mean；
4. `barycentric`；
5. `tree_local_factor`；
6. reliable-gene NNLS/ridge self-expression；
7. local sparse dictionary / LASSO；
8. graph harmonic 或 local linear；
9. local tree / random forest；
10. 若数据量足够，再比较轻量 local attention。

只有在 oracle masks、binomial thinning、DV preservation、marker leakage 和 uncertainty calibration 上同时稳定后，才应改变默认 recovery engine。

---

## 7. 关键参考文献

- Li WV, Li JJ. An accurate and robust imputation method scImpute for single-cell RNA-seq data. *Nat Commun*. 2018. doi:10.1038/s41467-018-03405-7.
- Chen M, Zhou X. VIPER: variability-preserving imputation for accurate gene expression recovery in single-cell RNA sequencing studies. *Genome Biol*. 2018. doi:10.1186/s13059-018-1575-1.
- Huang M et al. SAVER: Gene expression recovery for single-cell RNA sequencing. *Nat Methods*. 2018. doi:10.1038/s41592-018-0033-z.
- van Dijk D et al. Recovering gene interactions from single-cell data using data diffusion. *Cell*. 2018. doi:10.1016/j.cell.2018.05.061.
- Eraslan G et al. Single-cell RNA-seq denoising using a deep count autoencoder. *Nat Commun*. 2019. doi:10.1038/s41467-018-07931-2.
- Arisdakessian C et al. DeepImpute. *Genome Biol*. 2019. doi:10.1186/s13059-019-1837-6.
- Gong W et al. scTSSR. *Bioinformatics*. 2020. doi:10.1093/bioinformatics/btaa108.
- Wang J et al. scGNN. *Nat Commun*. 2021. doi:10.1038/s41467-021-22197-x.
- Linderman GC et al. Zero-preserving imputation of single-cell RNA-seq data / ALRA. *Nat Commun*. 2022. doi:10.1038/s41467-021-27729-z.
- Gong W et al. Imputing dropouts for single-cell RNA sequencing based on multi-objective optimization / scMOO. *Bioinformatics*. 2022. doi:10.1093/bioinformatics/btac300.
- Zhang H et al. scTsI. *Brief Bioinform*. 2025. doi:10.1093/bib/bbaf298.
- Huang S et al. D3Impute. *PLOS Comput Biol*. 2025. doi:10.1371/journal.pcbi.1013744.
- Mou Y, Li S, Ji G. scIRT. *Int J Mol Sci*. 2026. doi:10.3390/ijms27031173.
- Wu Y et al. Prior-guided factorization for reliable imputation of scRNA-seq data / scZN. *PLOS Comput Biol*. 2026. doi:10.1371/journal.pcbi.1014051.
- Vo LT et al. scZiva. *BMC Bioinformatics*. 2026. doi:10.1186/s12859-026-06422-2.
- scDDI. *Brief Bioinform*. 2026, bbag072.
- Borah K, Das HS. scGImpute. *Comput Biol Chem*. 2026. doi:10.1016/j.compbiolchem.2025.108856.
- Zi T et al. Stochastic hypergraph co-contrastive learning for single-cell RNA-seq data imputation. *Complex Intell Syst*. 2026. doi:10.1007/s40747-026-02360-x.
- Tian Y et al. GSCI: A Generative and Sparse Compressed Sensing Imputation Framework for Single-Cell RNA-Sequencing Dropout Recovery. *Comput Struct Biotechnol J*. 2026. doi:10.34133/csbj.0124.
- SCR-MF. arXiv:2511.16923, preprint.
- Iwashita Y et al. A Large-Scale Comparative Analysis of Imputation Methods for Single-Cell RNA Sequencing Data. arXiv:2603.24626, preprint.
