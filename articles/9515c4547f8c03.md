---
title: "被覆半径を効率よく小さくする点群の取り方について"
emoji: "☕"
type: "tech"  # tech or idea
topics: ["math", "機械学習", "最適化"]
published: true
---



# 導入

空間が与えられたとき，効率よく点群で近似することは，様々な場面で有用である．例えば，調べたい関数を点群の上で調べておいて，空間上の関数を近似するということが考えられる．今回は，空間を点群でどれくらい近似できているかを測定する被覆半径を紹介し，空間から一様に点群をとれば被覆半径を一番効率よく小さくできることを紹介する．


# 定義

距離空間$(X,\delta )$に対する部分集合$Y$の被覆半径とは

$$ \rho (X , Y) := \sup_{x \in X} \inf_{y \in Y} \delta (x,y) $$

である．被覆半径$\rho (X , Y)$は$X$と$Y$のHausdroff距離と見なすこともできる．

今回は，集合$X$から$n$個の有限点集合$X_n$をどのように取ると，被覆半径を効率よく小さくできるのかを考える．

以下，空間$(X,\delta)$は$C^1$級$d$次元コンパクトRiemann多様体とする．コンパクトRiemann多様体はEuclid空間にある有界でなめらかな図形，を想定している．

# 理論限界

[Reznikov-Saff][RS2016]は点群の被覆半径の理論限界を示した．
> 定理 1. \[[RS2016][RS2016], p. 3 \] 空間$(X,\delta)$は$C^1$級$d$次元コンパクトRiemann多様体であるとする．このとき，ある定数$C_1$が存在して，任意の正数$n \in \mathbb{Z}_{>0}$に対し，$n$点の任意の点群集合$X_n (\subset X)$は
> 
> $$ \rho (X , X_n ) \geq C_1 n^{-1/d}  .$$
> 

※実際はより一般の場合で示している．

証明の概要である．測度$\mu$をRiemann計量から誘導される測度とし．$B(x,r)$を点$x \in X$の半径$r(>0)$の球とする．実数$r >0 $は

$$ \mu (B(x,r) ) \leq C (2r)^d \tag{1} $$

を満たす．ここで，$X_n = \\{ y_1 , \ldots , y_n \\}$として，

$$ \rho_n = \rho ( X , X_n )$$

とする．このとき，被覆半径$\rho_n$の定義と(1)から

$$ 0 < \mu (X) \leq \sum_{i=1}^n \mu (B(y_i , \rho_n ) ) \leq n C (2 \rho_n )^d $$

であるから，定理が従う．


# 一様に点群を取る


[Flatto-Newman][FN1977]によって，定理1の理論限界を達成する以下の定理が示された．
> 定理 2. \[[FN1977][FN1977], 定理 2.1\]
> 空間$(X,\delta)$は$C^1$級$d$次元コンパクトRiemann多様体であるとする．このとき，ある定数$C_2$が存在して，任意の正数$n \in \mathbb{Z}_{>0}$に対し，$n$点の点群集合$X_n (\subset X)$で以下を満たすものが存在する：
> 
> $$ \rho (X , X_n ) \leq C_2 n^{-1/d}  .$$
> 

証明の概要についてである．まず簡単のために，超立方体$X=[ 0 , 1 ]^d$上のEuclid距離の場合で考える．
このとき，$X_{n^d}$を以下で取る．

$$ X_{n^d} = \left\{ \left( \frac{i_1}{n-1} , \ldots , \frac{i_d}{n-1} \right) \ \middle| \ i_1 , \ldots i_d = 0 , \ldots n-1 \right\}  $$

これは，

$$  \rho ([0,1]^d , X_{n^d} ) \leq C_2 n $$

であることが従う．Euclid距離でない場合も，適切に距離を評価することによって，得ることができる．

一般の場合については，$X$のアトラス$\{ (U_i = [0,1]^d ,\phi_i )\}_{i \in I}$であって，$X$はコンパクトより$I$は有限集合で取れる．各局所座標系$(U_i = [0,1]^d ,\phi_i )$で点群をとり，$I$について考えることで，欲しい物を得る．

# 関数の点群近似について

空間$(X,\delta)$からノルム空間への写像$f$を点群$X_n$を用いて近似することを考えた際の近似精度は，[Ferber-Huang-Zha達][FHZ+2023]によって研究されてた．

> 定理 3. \[[FHZ+2023][FHZ+2023], 補題 4.2\]
> 写像$f$は$L$-Lipschitzとする，すなわち，任意の$x, y \in X$に対して
>
> $$ \| f (x ) - f(y) \| \leq L \delta (x,y) $$
> を満たすとする．点群$X_n$で近似した写像$\check{f}\_n $を
>
> $$ \check{f}\_n (x) := f (\check{x} ) , \ \text{ここで} \check{x} \in \mathrm{arg \ min}\_{y \in X\_n } \delta (x , y) $$
>
> と定義する．このとき，$f$と$\check{f}\_n$の差は
> 
> $$\sup_{x \in X} \| f (x) - \check{f}_n (x) \| \leq L C_2 n^{-1/d}$$
> 

証明

$$
\begin{align*}
    \sup_{x \in X} \| f (x) - \check{f}_n (x) \|
    & = \sup_{x \in X} \| f (x) - f (\check{x}) \| \\
    & \leq
    \sup_{x \in X} L \delta (x ,\check{x}) \\
    & =
    \sup_{x \in X} L \inf_{y \in X_n} \delta (x ,y)
    \\
    & =
    L \rho (X  ,X_n )
    \leq L C_2 n^{-1/d}
\end{align*}
$$

となる．

# 結論

定理1, 2から，空間から一様に点群をとれば被覆半径を一番効率よく小さくできる，$O(n^{-1/d})$で収束する．

また，空間から一様に点群を取って，点群で関数を近似する際の精度は上から，$O(n^{-1/d})$で抑えられる．

# 参考文献

\[[FN1977][FN1977]\] L. Flatto, D. J. Newman, Random coverings, Acta Math. 138: 241-264 (1977). DOI: 10.1007/BF02392317

\[[RS2016][RS2016]\] A. Reznikov, E. B. Saff, The Covering Radius of Randomly Distributed Points on a Manifold, International Mathematics Research Notices, Volume 2016, Issue 19, 2016, Pages 6065–6094, https://doi.org/10.1093/imrn/rnv342

\[[FHZ+2023][FHZ+2023]\] Aaron M Ferber, Taoan Huang, Daochen Zha, Martin Schubert, Benoit Steiner, Bistra Dilkina, Yuandong Tian, SurCo: Learning Linear SURrogates for COmbinatorial Nonlinear Optimization Problems, Proceedings of the 40th International Conference on Machine Learning, PMLR 202:10034-10052, 2023.

[FN1977]:https://projecteuclid.org/journals/acta-mathematica/volume-138/issue-none/Random-coverings/10.1007/BF02392317.full

[RS2016]:https://academic.oup.com/imrn/article-abstract/2016/19/6065/2399425?redirectedFrom=fulltext&login=true

[FHZ+2023]:https://proceedings.mlr.press/v202/ferber23a.html
