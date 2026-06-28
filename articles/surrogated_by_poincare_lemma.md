---
title: "サロゲートモデル作成のアイデア ～Poincaréの補題を用いて～"
emoji: "🔑"
type: "idea"  # tech or idea
topics: ["math", "機械学習", "微分幾何", "微積分", "数理最適化"]
published: true
---


# この記事を読んで得られること
- ベクトル値関数のパラメータを推定するために，Poincaréの補題（微積分の基本定理）を用いると，サロゲート関数を構成できる場合があることがわかる．
- 例として，逆最適化問題の可解性を判定するときに，このアイデアを用いると，凸なサロゲート関数（suboptimality損失）を構成でき，最小値を求める問題に帰着できる話をする．

　
# 逆最適化問題の定義

（整数線形計画や混合整数線形計画に現れる）実行可能領域を$X \subset \mathbb{R}^d$とする．重みを$\theta \in \Delta^{d-1} = \{ \theta \in \mathbb{R}_{\geq 0}^d \mid \sum_i \theta_i = 1 \}$とする．順問題を

$$
    x^* (\theta) \in \argmin_{x \in X} \theta^\top x = \sum_{i=1}^d \theta_i x_i
$$

とする．すなわち，重み$\theta$が与えられるので最適解$x^* (\theta)$を求めよという問題である．このとき，逆最適化問題を，未知の重み$\theta^* \in \Delta^{d-1}$を用いて$\hat{x} = x^{*} (\theta^*)$と表されるとしたとき，この未知の重み$\theta^*$を求める問題とする，すなわち，

$$
    \hat{x} \in \argmin_{x \in X} \theta^\top x = \sum_{i=1}^d \theta_i x_i
$$

となる$\theta \in \Delta^{d-1}$を求める問題である．


# 逆最適化問題を予測損失最小化で解こうとしたときの課題

逆最適化問題のナイーブな定式化として，予測誤差最小化問題として定式化することが挙げられる，すなわち，

$$
    \text{minimize } \ell_{\mathrm{pre}} (\theta) := \| x^* (\theta) - \hat{x} \|
    \text{ subject to } \theta \in \Delta^{d-1}
$$

と定式化できる．しかし，（整数線形計画であらわれる）実行可能領域の場合，$x^* (\theta) \in X$は離散的な集合$X$の値しかとらないので，勾配法を使って求めるということが難しい．


# Poincaréの補題

Poincaréの補題とは，ベクトル値関数$g \colon \mathbb{R}^d \to \mathbb{R}^d$が与えられたとき，$\nabla f = g$となる関数$f \colon \mathbb{R}^d \to \mathbb{R}$が存在するかを判定し，存在すれば具体的に構成できるものである．

> **補題 1.**　(Poincaréの補題 (cf. \[Tsuboi 2008\]))　ベクトル値関数$g = (g_1 ,\ldots , g_d) \in C^1 ( \mathbb{R}^d , \mathbb{R}^d)$とする．このとき，$\nabla f = g$となる関数$f \in C^1 ( \mathbb{R}^d , \mathbb{R})$が存在する場合，$f$は
>
> $$
    f(\theta) - f(0) = \int_0^1 \sum_{i=1}^d \theta_i g_i ( s \theta ) ds
$$
> 
> で与えられる．


**証明**

$$
    f (\theta) - f(0) = \int_0^1 \frac{d f ( s \theta)}{d s} ds 
$$

となる．微分の連鎖律より，

$$
    \int_0^1 \frac{d f ( s \theta)}{d s} ds
    =
    \int_0^1 \sum_{i=1}^d \theta_i \frac{\partial f ( s \theta)}{\partial \theta_i} ds
    = 
    \int_0^1 \sum_{i=1}^d \theta_i g_i ( s \theta) ds
$$

となる．まとめると，

$$
    f (\theta) - f (0) = \int_0^1 \sum_{i=1}^d \theta_i g_i ( s \theta) ds 
$$

となる．Q.E.D.



# Poincaréの補題とサロゲート関数

Poincaréの補題の仮定を満たす$f$が存在して，さらに，$f$がLipschitz連続な凸関数で，最小値が存在するとする．
このとき，以下の問題は同値であると考えられる．

$$
    g(\theta) = 0
    \quad 
    \Leftrightarrow
    \quad 
    \text{minimize } \| g(\theta) \|
    \quad 
    \Leftrightarrow
    \quad 
    \text{minimize } f (\theta)
$$

この考察は$g(\theta) = 0$を解く問題を，サロゲートして，$\text{minimize } f (\theta)$を解く問題に帰着させるアイデアである．


# Poincaréの補題のアイデアを逆最適化問題に適用する

逆最適化において，$g (\theta ) = x^* (\theta) - \hat{x}$，もしくは，逆最適化の予測損失$\ell_{\mathrm{pre}} = \| g \|$として，Poincaréの補題のアイデアを用いると，サロゲート関数$\ell_{\mathrm{sub}}$は

$$
\begin{align*}
    \ell_{\mathrm{sub}} (\theta) 
    & = \int_0^1 \sum_{i=1}^d \theta_i (x^* ( s \theta) - \hat{x}) ds \\
    & = \int_0^1 \sum_{i=1}^d \theta_i (x^* ( \theta) - \hat{x}) ds \\
    & = \int_0^1 \theta^\top (x^* (\theta) - \hat{x}) ds \\
    & = \theta^\top (x^* (\theta) - \hat{x})
\end{align*}
$$

となる．上記のサロゲート関数をsuboptimality損失と呼ぶ．suboptimality損失はLipschitz連続な凸関数であり，その（劣）勾配が$x^* (\theta) - \hat{x}$になることが知られている（[Barman+ 2018][Barman2018], 命題 3.1）．このことから，以下は（ほぼ）同値だと考えられる．

$$
    x^* (\theta) = \hat{x}
    \quad 
    \Leftrightarrow
    \quad 
    \text{minimize } \ell_{\mathrm{pre}} (\theta)
    \quad 
    \Leftrightarrow
    \quad 
    \text{minimize } \ell_{\mathrm{sub}} (\theta)
$$

このことから，suboptimality損失を最小化すれば，逆最適化問題は解けると考えられる．Lipschitz連続な凸関数を最小値に近づけるアルゴリズムとして，射影劣勾配法が知られている．実は，以下が成り立つ，

> **定理 3.**　(cf. [Kitaoka 2024][Kitaoka2024], 定理 4.11, 4.12)　最適解$x^* (\theta^*)$は一意に定まるとする．$\{ \theta^t \}$をsuboptimality損失$\ell_{\mathrm{sub}}$に関して学習率を$\alpha t^{-1/2}$もしくはステップ幅を$\alpha t^{-1/2}$とする射影劣勾配法で得られた点列とする．このとき，十分大きい$t \gg 1$で，$x^* (\theta^t) = \hat{x}$，つまり，逆最適化問題は完全に解ける．



# 参考文献

[Tsuboi 2008] 坪井俊, 幾何学 III 微分形式 大学数学入門 6, 東京大学出版会, (2008). 

\[[Barman+ 2018][Barman2018]\] Andreas Bärmann, Alexander Martin, Sebastian Pokutta, Oskar Schneider, An Online-Learning Approach to Inverse Optimization, arXiv:1810.12997. 

\[[Kitaoka 2024][Kitaoka2024]\] Akira Kitaoka, Exact Solution to Data-Driven Inverse Optimization of MILPs in Finite Time via Gradient-Based Methods, arXiv:2405.14273v7. 

[Barman2018]:https://arxiv.org/abs/1810.12997

[Kitaoka2024]:https://arxiv.org/abs/2405.14273v7
