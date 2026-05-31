---
title: "ランダムサンプリングによる被覆半径の期待値"
emoji: "🎲"
type: "tech"  # tech or idea
topics: ["math", "機械学習", "確率", "確率論"]
published: true
---

# 導入

空間が与えられたとき，ランダムに点を取ることにより，どれくらい空間を近似できるかを考える．今回は，空間からランダムに点を取ったときにどれくらい近似できているかを測定する被覆半径を紹介し，被覆半径の値の確率とその期待値について紹介する．

# 被覆半径の定義

距離空間$( X,\delta )$に対する部分集合$Y$の被覆半径とは

$$ \rho (X , Y) := \sup_{x \in X} \inf_{y \in Y} \delta (x,y) $$

である．被覆半径$\rho (X , Y)$は$X$と$Y$のHausdroff距離と見なすこともできる．

今回は，集合$X$から$n$個ランダムに取ったときの，被覆半径は何なのか，特に，期待値はいくらなのかについて考える．

以下，空間$(X,\delta)$は$C^1$級$d$次元コンパクト(角付き)Riemann多様体とする．コンパクトRiemann多様体はEuclid空間にある有界でなめらかな図形，を想定している．

# ランダムサンプリング

[Reznikov-Saff][RS2016]によって，被覆半径の確率評価，及び，期待値について，計算されている．

> 定理 1. \[[RS2016][RS2016] , 系2.3 \] 空間$(X,\delta)$は$C^1$級$d$次元コンパクト(角付き)Riemann多様体とする．また，$X_n$を空間$(X,\delta)$上の体積測度由来の確率分布から$n$点サンプリングした確率分布とする．このとき，ある$c_1, c_2 >0$が存在して，任意の$\varepsilon >0$に対して，ある$N(\varepsilon) \in \mathbb{Z}_{>0}$が存在して，$N > N(\varepsilon)$ならば
> 
> $$
    \mathbb{P} \left[ c_1 \left( \frac{\log N}{N} \right)^{1/d} \leq \rho (X_N , X) \leq c_2 \left( \frac{\log N}{N} \right)^{1/d} \right] > 1 - \varepsilon .$$
> 
> その上，ある定数$C_1, C_2 > 0$が存在して，
> 
> $$ C_1 \left( \frac{\log N}{N} \right)^{1/d} \leq \mathbb{E} \rho (X_N , X) \leq C_2 \left( \frac{\log N}{N} \right)^{1/d} . $$
> 

※実際はより一般の場合で示している．

証明は[RS2016][RS2016]に投げることとする．

# 理論値限界との比較

[Reznikov-Saff][RS2016]は点群の被覆半径の理論限界を示した．
> 定理 2. \[[RS2016][RS2016], p. 3 \] 空間$(X,\delta)$は$C^1$級$d$次元コンパクトRiemann多様体であるとする．このとき，ある定数$C_1$が存在して，任意の正数$n \in \mathbb{Z}_{>0}$に対し，$n$点の任意の点群集合$X_n (\subset X)$は
> 
> $$ \rho (X , X_n ) \geq C_1 n^{-1/d} . $$
> 

※実際はより一般の場合で示している．

定理 1.より，ランダムに取った場合より，$\log$の分だけ計算効率が悪いことが言える．
また，[Flatto-Newman][FN1977]により，定理2の不等号を定数倍を除いて達成する点列の取り方が存在する[FN1977][FN1977]．つまり，最良に点を取ったときに比べて，$\log$の分だけ計算効率が悪いことが言える．

# 参考文献

\[[FN1977][FN1977]\] L. Flatto, D. J. Newman, Random coverings, Acta Math. 138: 241-264 (1977). DOI: 10.1007/BF02392317

\[[RS2016][RS2016]\] A. Reznikov, E. B. Saff, The Covering Radius of Randomly Distributed Points on a Manifold, International Mathematics Research Notices, Volume 2016, Issue 19, 2016, Pages 6065–6094, https://doi.org/10.1093/imrn/rnv342

[FN1977]:https://projecteuclid.org/journals/acta-mathematica/volume-138/issue-none/Random-coverings/10.1007/BF02392317.full

[RS2016]:https://academic.oup.com/imrn/article-abstract/2016/19/6065/2399425?redirectedFrom=fulltext&login=true
