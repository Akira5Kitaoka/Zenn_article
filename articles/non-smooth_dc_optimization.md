---
title: "滑らかでないDC最適化"
emoji: "📉"
type: "idea"  # tech or idea
topics: ["数理最適化", "DC最適化"]
published: false
---

# これは何？
空間$X \subset \mathbb{R}^d$を有界凸集合とする．関数$g, \, h \colon X \to \mathbb{R}$を（滑らかとは限らない）Lipschitz関数とする．滑らかでないDC関数 $f := g-h$ (difference of convex functions)としたとき，以下の最適化問題をDC最適化問題という．

$$
    \min_{x \in X} f(x) = \min_{x \in X} \left( g(x) - h(x) \right)
$$

今回はDC最適化問題を解く方法，DC最適化手法について解説する．DC関数$f$は凸関数でないため，$f$の最小値を得ることは難しいと考えられている．その代わり，$f$の臨界点に漸近するアルゴリズムとして，切除平面の選び方に近接バンドル法を活用したDC最適化アルゴリズム\[[dO2019][dO2019]\]について紹介する．

# 問題設定
- $X$は有界凸集合．
- 関数$g, \, h \colon X \to \mathbb{R}$のオラクルは計算できる．
- (劣)微分(の元を与えるもの)$\nabla g, \, \nabla h \colon X \to \mathbb{R}^d$のオラクルは計算できる．

# DCアルゴリズム（DCA）の型

## アルゴリズム

1. 初期値$x^0 \in X$
1. for $k = 0 , \ldots , K-1$
1. $\quad$ $x^{k+1} \gets \argmin_{x \in X} (g(x) - \langle \nabla h(x^{k}) , x \rangle)$
1. endfor
1. return $x \in \{ x^1 , \ldots , x^K \}$ such that $f(x) = \min_{k = 1 , \ldots , K} f(x^k)$

各反復で劣勾配$\nabla h(x^{k})$を使って凸最適化問題を解く，というアルゴリズムである．

## DCAの課題
- 各反復において，凸最適化問題をどれくらい近似的に解くべきか考えないといけない．一般に凸最適化問題を厳密に解くことは難しい．
- 有界凸集合$X$と，オラクル$g, h, \nabla g, \nabla h$のみを使える形になっていない．

# de Oliveira 2019：DC最適化のための近接バンドル法

\[[dO2019][dO2019]\]は，DCAの課題のうち1つ目——毎反復の凸最適化問題を厳密に解かねばならない——を，**毎反復ちょうど1個の凸2次計画に置き換える**手法である．$g$を切除平面モデル（区分線形な下界）で近似し，安定中心のまわりに解を引き留める近接項を加える．

## 記号の準備

**カット．** 点$x^j$でオラクルを呼ぶと組$(g(x^j), \nabla g(x^j))$が手に入る．$g$の凸性から，この組はアフィン関数による$g$の下界

$$
    \bar{g}^j(x) := g(x^j) + \langle \nabla g(x^j) , x - x^j \rangle \le g(x) \quad (\forall x \in X)
$$

を1本与える．これを**カット**と呼ぶ．添字$j \ge 0$は点$x^j$に対応する．

**近接項．** 2回微分可能な強凸関数$\omega$を1つ選び，Bregman関数$D(x,y) := \omega(x) - \omega(y) - \langle \nabla \omega(y) , x-y \rangle$を近接項に使う．$\omega$は$g, h$とは無関係に**こちらが自由に選ぶ**関数なので，$g, h$がLipschitzなだけでも問題ない．以下では標準的な選択$\omega(\cdot) = \| \cdot \|_2^2/2$，すなわち

$$
    D(x,y) = \frac{1}{2} \| x - y \|_2^2
$$

を採用する（論文の数値実験もこの選択である）．

**降下判定と安定中心．** アルゴリズムは反復点$x^0, x^1, x^2, \ldots$を生成するが，そのすべてを近接項の中心に使うわけではない．中心に採用してよいのは，目的関数を十分に減らした点だけである．この「十分に減らした」を判定する不等式が次の$(D)$であり，これによって安定中心が定まる．なお$\kappa \in (0,1)$と$\underline{\mu} > 0$はアルゴリズムの入力パラメータである．

> 定義 1. \[[dO2019][dO2019], 式 (13)\] 反復$k$の開始時点での安定中心の番号を$\ell_k$と書く．$k(0) := 0$，$\ell_0 := 0$とし，反復$k$で得られた点$x^{k+1}$に対し
>
> $$ g(x^{k+1}) \le g(x^{k(\ell_k)}) + \langle \nabla h(x^{k(\ell_k)}) , x^{k+1} - x^{k(\ell_k)} \rangle - \frac{\kappa \underline{\mu}}{2} \| x^{k+1} - x^{k(\ell_k)} \|_2^2 \tag{D}$$
>
> が成り立つ場合を**serious step**と呼び，$\ell_{k+1} := \ell_k + 1$，$k(\ell_{k+1}) := k+1$と定める．$(D)$が成り立たない場合を**null step**と呼び，$\ell_{k+1} := \ell_k$と定める（$k(\cdot)$は更新しない）．こうして定まる点$x^{k(\ell)}$を第$\ell$**安定中心(stability center)**と呼ぶ．

$k(\ell)$は第$\ell$安定中心が得られた反復番号であり，$\ell_k$は反復$k$の時点での最新の安定中心の番号である．混乱の恐れがないときは$\ell_k$を単に$\ell$と書き，反復$k$における安定中心を$x^{k(\ell)}$と書く．$(D)$は$f = g-h$ではなく$g$だけで書かれている点に注意する（この理由は後述する）．

## 部分問題

反復$k$の開始時点で，アルゴリズムは次の2つを持っているものとする．

- 有限集合$B^k$（**バンドル**と呼ぶ）
- 各$j \in B^k$に対する，$g$のアフィンな下界$\bar{g}^j \le g$

出発点は$B^0 := \{ 0 \}$であり，$\bar{g}^0$は上で定義したカットである．反復が進んでもこの形が保たれることは，$k$についての帰納法で従う（次節で確認する）．バンドルから，$g$の区分線形な下界

$$
    \check{g}^k(x) := \max_{j \in B^k} \bar{g}^j(x) \le g(x)
$$

が定まる．これを**切除平面モデル**と呼び，部分問題の中で$g$の代わりに使う．$|B^k|$が大きいほどモデルは$g$に近づくが，部分問題は重くなる．

安定中心$x^{k(\ell)}$と近接パラメータ$\mu_k > 0$に対し，次の狭義凸計画を解く．

$$
    \min_{x \in X, \, r \in \mathbb{R}} \ r - \langle \nabla h(x^{k(\ell)}) , x - x^{k(\ell)} \rangle + \frac{\mu_k}{2} \| x - x^{k(\ell)} \|_2^2
    \quad \text{s.t.} \quad \bar{g}^j(x) \le r , \ \ j \in B^k
    \tag{P}
$$

補助変数$r$と$|B^k|$本の1次不等式は，$\max$で書かれた切除平面モデル$\check{g}^k$を線形制約に開いたものであり，最適解では$r = \check{g}^k(x)$になる．DCAと同じく$h$は安定中心での線形化$\langle \nabla h(x^{k(\ell)}), \cdot \rangle$で置き換えているが，$g$のほうも切除平面モデルに置き換わっている点が違う．$X$が多面体なら$(P)$は凸2次計画（QP）である．

$(P)$は$x$について狭義凸なので最適解が一意に定まる．これを$(x^{k+1}, r^{k+1})$と書き，制約$\bar{g}^j(x) \le r$に対するLagrange乗数を$\alpha_j \ge 0$（$j \in B^k$）と書く（乗数の存在には$X$が多面体であるか，適当な制約想定が成り立つことを仮定する）．Lagrange関数を$r$で偏微分すると$1 - \sum_{j \in B^k} \alpha_j = 0$，すなわち

$$
    \sum_{j \in B^k} \alpha_j = 1
$$

が最適性条件から従う．つまり$(\alpha_j)_{j \in B^k}$は凸結合の係数である．

## 集約カット

$(P)$を1回解くたびに，バンドル全体を1本に圧縮したカットが手に入る．上の乗数による凸結合

$$
    \bar{g}^{-k}(x) := \sum_{j \in B^k} \alpha_j \, \bar{g}^j(x)
$$

を考えると，アフィン関数の凸結合なので$\bar{g}^{-k}$はアフィンであり，その勾配は$p^{k+1} := \sum_{j \in B^k} \alpha_j \nabla g(x^j)$である．また各$\bar{g}^j$が$g$の下界であることと$\alpha_j \ge 0$，$\sum_j \alpha_j = 1$から

$$
    \bar{g}^{-k}(x) = \sum_{j \in B^k} \alpha_j \, \bar{g}^j(x) \le \sum_{j \in B^k} \alpha_j \, g(x) = g(x)
$$

となるので，$\bar{g}^{-k}$も$g$のアフィンな下界である．これを**集約カット**と呼び，添字$-k$を割り当てる．

:::message
相補性条件より，$\alpha_j > 0$となる$j$では$\bar{g}^j(x^{k+1}) = r^{k+1} = \check{g}^k(x^{k+1})$である．したがって$\bar{g}^{-k}(x^{k+1}) = \check{g}^k(x^{k+1})$となり，$\bar{g}^{-k}$は切除平面モデル$\check{g}^k$に点$x^{k+1}$で接する：

$$
    \bar{g}^{-k}(x) = \check{g}^k(x^{k+1}) + \langle p^{k+1} , x - x^{k+1} \rangle
$$

論文ではこちらの形で書かれている．$g$の下界であることは$p^{k+1} \in \partial \check{g}^k(x^{k+1})$と$\check{g}^k \le g$からも言えるが，上の凸結合による議論のほうが短い．
:::

**利用できるカット．** 反復$k$の終わりに手元にあるカットの添字は

$$
    \{ 0, 1, \ldots, k+1 \} \cup \{ -k \}
$$

である．どの添字$j$についても$\bar{g}^j$がこの時点で確定しているからである：$j \le k$のカットは過去のオラクル呼び出しから，$j = k+1$のカットは反復$k$で計算した$(g(x^{k+1}), \nabla g(x^{k+1}))$から，$j = -k$は上の集約カットから得られる．次の反復で使うバンドル$B^{k+1}$は，この集合の部分集合として選ぶ．

## アルゴリズム

\[[dO2019][dO2019], Algorithm 1\]である（以下**アルゴリズム1**と呼ぶ）．

**入力**：$x^0 \in X$，$\kappa \in (0,1)$，$0 < \underline{\mu} \le \mu_0 \le \bar{\mu} < \infty$，許容誤差$\delta_{\mathrm{Tol}} \ge 0$

1. $g(x^0), \, \nabla g(x^0), \, \nabla h(x^0)$を計算し，$B^0 \gets \{ 0 \}$，$k(0) \gets \ell \gets 0$
1. for $k = 0, 1, 2, \ldots$
1. $\quad$ $(P)$を解いて$x^{k+1}$と乗数$\alpha_j \ge 0$を得る
1. $\quad$ **［停止判定］** if $\| x^{k+1} - x^{k(\ell)} \|_2 \le \delta_{\mathrm{Tol}}$ then return $x^{k(\ell)}$
1. $\quad$ **［オラクル］** $g(x^{k+1})$と$\nabla g(x^{k+1})$を計算
1. $\quad$ if $(D)$ が成り立つ
1. $\quad$ $\quad$ **［serious step］** $\nabla h(x^{k+1})$を計算し，定義 1. に従って$x^{k(\ell+1)} \gets x^{k+1}$，$k(\ell+1) \gets k+1$，$\ell \gets \ell + 1$
1. $\quad$ $\quad$ $B^{k+1}$を $\{ k+1 \} \subset B^{k+1} \subset \{ 0, 1, \ldots, k+1 \} \cup \{ -k \}$ を満たすように選び，$\mu_{k+1}$を区間$[ \underline{\mu} , \mu_k ]$から選ぶ
1. $\quad$ else
1. $\quad$ $\quad$ **［null step］** 定義 1. に従い，安定中心と$\ell$は据え置く
1. $\quad$ $\quad$ $B^{k+1}$を $\{ -k , \, k(\ell) , \, k+1 \} \subset B^{k+1} \subset \{ 0, 1, \ldots, k+1 \} \cup \{ -k \}$ を満たすように選び，$\mu_{k+1}$を区間$[ \mu_k , \bar{\mu} ]$から選ぶ
1. $\quad$ endif
1. endfor

降下判定$(D)$は，$f = g-h$そのものではなく$g$だけで書かれている．$h$の凸性から，これが成り立てば

$$
    f(x^{k(\ell+1)}) \le f(x^{k(\ell)}) - \frac{\kappa \underline{\mu}}{2} \| x^{k(\ell+1)} - x^{k(\ell)} \|_2^2
$$

が従うので，安定中心の列に沿って$f$は確かに減少する．

### バンドルと$\mu_k$の選び方の自由度

8行目・11行目は$B^{k+1}$と$\mu_{k+1}$を**一意に決めていない**．満たすべき条件だけを課している．つまりこの擬似コードは1本のアルゴリズムではなく**アルゴリズムの族**を定めており，条件を守る限りどう選んでも収束定理が成り立つ．実装者はメモリ量や部分問題の重さに応じて自由に決められる．

**バンドルの両端．** 上限$\{ 0, 1, \ldots, k+1 \} \cup \{ -k \}$は「使えるカットの全体」，下限は「最低限これは残せ」である．両端がそのまま2つの極端な実装に対応する．

- **最小**：$B^{k+1} \gets \{ -k , k(\ell) , k+1 \}$（null step）／$B^{k+1} \gets \{ k+1 \}$（serious step）．バンドルは常に3本・1本で，部分問題は最小サイズに保たれる．
- **最大**：$B^{k+1} \gets \{ 0, 1, \ldots, k+1 \} \cup \{ -k \}$．カットを一切捨てない．モデルは最も精密だが，部分問題が反復ごとに重くなる．

実用上はこの中間を取り，バンドルの上限本数を決めて古いカットや非活性なカットを捨てる．

**null stepで3本が必要な理由．** 3本はそれぞれ役割が違い，どれか1つでは代用できない．

| 添字 | カット | 役割 |
| --- | --- | --- |
| $k+1$ | 試行点$x^{k+1}$で計算した新しいカット | $x^{k+1}$の近くでモデル$\check{g}$を改善する．これがないと次の反復も同じ$x^{k+1}$が返ってきて進まない |
| $k(\ell)$ | 現在の安定中心$x^{k(\ell)}$のカット | $\check{g}^k(x^{k(\ell)}) = g(x^{k(\ell)})$を保証する．停止判定が臨界性を意味することの証明に直接使う |
| $-k$ | 集約カット | 捨てた過去のカットの情報を1本に圧縮して保持する |

null stepは「モデルが粗すぎて良い点が出なかった」場合なので，モデルを改善しつつ過去の情報を失わない，という要請である．一方serious stepは「モデルが十分に良かった」場合なので，バンドルをリセットしてよい．$\{ k+1 \}$は，serious stepでは$k(\ell+1) = k+1$すなわち$x^{k+1}$が新しい安定中心になることに注意すれば，上の表の2行目（新しい安定中心のカット）と同じものである．

**$\mu_k$の両端．** $\mu_k$の更新規則は非対称で，serious stepでは下げてよく（$\mu_{k+1} \le \mu_k$），null stepでは上げるだけ（$\mu_{k+1} \ge \mu_k$）である．うまく進めたら近接項を弱めて大きく動けるようにし，失敗したら近接項を強めて安定中心の近くに引き留める，という意図である．例えばserious stepで$\mu_{k+1} = \max \{ \underline{\mu} , \mu_k / 2 \}$，null stepで$\mu_{k+1} = \min \{ \bar{\mu} , 2 \mu_k \}$とすればよい．

両端$\underline{\mu}, \bar{\mu}$は飾りではない．

- $\underline{\mu} > 0$がないと$\mu_k \to 0$が許される．降下判定が保証する1回あたりの減少量は$\frac{\kappa \underline{\mu}}{2} \| x^{k+1} - x^{k(\ell)} \|_2^2$なので，$\underline{\mu} = 0$では$f$の減少の保証が消える．
- $\bar{\mu} < \infty$がないと$\mu_k \to \infty$が許される．近接項が強すぎると$x^{k+1} \to x^{k(\ell)}$となるので，臨界点に近づいていないのに停止判定$\| x^{k+1} - x^{k(\ell)} \|_2 \le \delta_{\mathrm{Tol}}$が発火してしまう．

**null step後に$\mu_k$を減らすことが禁止**されている点も重要で，これが収束証明の要になっている．

## 理論保証

まず，収束先を述べるための概念を用意する．$\bar{x} \in X$における$X$の**法錐**を

$$
    N_X(\bar{x}) := \{ s \in \mathbb{R}^d : \langle s , x - \bar{x} \rangle \le 0 \quad \forall x \in X \}
$$

とする．

> 定義 2. \[[dO2019][dO2019], 式 (3)\] 点$\bar{x} \in X$が$f(x)$の$X$に関する**臨界点(critical point)**であるとは，
>
> $$ \partial h(\bar{x}) \cap \left( \partial g(\bar{x}) + N_X(\bar{x}) \right) \neq \emptyset \tag{C} $$
>
> が成り立つことをいう．

DC最適化では最適性の概念に強さの順序があり，

$$
    \text{局所最小解} \Rightarrow \text{臨界点}
$$

が成り立つ（\[[dO2019][dO2019]\] 2節）．逆向きは一般には成り立たない．臨界点は弱い概念だが，$X$と4つのオラクルだけで検証しうる形をしている点が実用上の利点である．

以下，$\mathcal{L} \subset \{ 0, 1, 2, \ldots \}$をserious stepの番号の集合，すなわち$\ell \in \mathcal{L}$のとき$x^{k(\ell)}$が$\ell$番目の安定中心であるとする．また$i_X$を$X$の指示関数（$x \in X$で$0$，そうでなければ$+\infty$），$\partial_\epsilon \phi(x) := \{ s : \phi(y) \ge \phi(x) + \langle s , y - x \rangle - \epsilon \ \ \forall y \}$を$\epsilon$-劣微分，$B(0;\rho)$を原点中心・半径$\rho$の閉球とする．

> 定理 3. \[[dO2019][dO2019], Theorem 1\] 次を仮定する：$X \subset \mathbb{R}^d$は空でない閉凸集合で，ある開集合$\Omega \supset X$上で$g, h$は有限値の凸関数である；準位集合$\{ x \in X : f(x) \le f(x^0) \}$は有界である；$\omega$は$X$上で2回微分可能かつ強凸であり，すべての$k$で$0 < \underline{\mu} \le \mu_k \le \bar{\mu} < \infty$である；上のアルゴリズムのバンドルと$\mu_k$の更新規則が守られている．このとき次が成り立つ．
>
> (i) $\delta_{\mathrm{Tol}} = 0$ならば，生成される安定中心の列$\{ x^{k(\ell)} \}_{\ell \in \mathcal{L}}$の任意の集積点$\bar{x}$は$(C)$を満たす，すなわち臨界点である．
>
> (ii) $\delta_{\mathrm{Tol}} > 0$ならば，アルゴリズムは有限回の反復で停止する．さらに$\nabla \omega$が局所Lipschitz連続ならば，定数$s_1, s_2 > 0$が存在して，出力$\bar{x} = x^{k(\ell)}$は
>
> $$ \partial_{s_1 \delta_{\mathrm{Tol}}} \left[ g(\bar{x}) + i_X(\bar{x}) \right] \cap \left[ \partial h(\bar{x}) + B(0 ; s_2 \delta_{\mathrm{Tol}}) \right] \neq \emptyset $$
>
> を満たす（近似臨界点）．

凸性から$\partial [ g(\bar{x}) + i_X(\bar{x}) ] = \partial g(\bar{x}) + N_X(\bar{x})$が成り立つので，(ii)の包含式で形式的に$\delta_{\mathrm{Tol}} = 0$とおくと$(C)$そのものになる．つまり(ii)は$(C)$の$\delta_{\mathrm{Tol}}$-緩和版である．定数は$s_1 = M + L$（$M, L$は劣勾配の有界性から来る未知定数），$s_2 = \bar{\mu} L_\omega$（$L_\omega$は$\nabla \omega$のLipschitz定数）であり，本記事の選択$\omega = \| \cdot \|_2^2 / 2$では$\nabla \omega = \mathrm{id}$なので$L_\omega = 1$，したがって$s_2 = \bar{\mu}$である．

<!--
:::message
注釈
定理 3. の(i)は，反復の生成のされ方による3つの場合分けからなり，それぞれ別の補題・命題が担当する．

| 場合 | 担当 |
| --- | --- |
| 有限回で停止する | \[[dO2019][dO2019], Lemma 2\]：停止時は$x^{k+1} = x^{k(\ell)}$である．$k(\ell) \in B^k$より$\check{g}^k(x^{k(\ell)}) = g(x^{k(\ell)})$なので，この点はモデル$\check{g}^k$を$g$に戻した真の部分問題も解いており，その最適性条件が$(C)$を与える |
| 停止せず，serious stepが無限回 | \[[dO2019][dO2019], Proposition 2\]：$f$の減少量の総和が有限であることから$\| x^{k(\ell+1)} - x^{k(\ell)} \|_2 \to 0$を出し，$\epsilon$-劣微分の極限を取る |
| 停止せず，serious stepが有限回（以降null stepのみ） | \[[dO2019][dO2019], Proposition 3\]：最後の安定中心$\hat{x}$以降は$\nabla h(\hat{x})$が固定されるため通常の凸バンドル法に帰着し，$x^{k+1} \to \hat{x}$かつ$\hat{x}$が臨界点 |

3番目の場合の\[[dO2019][dO2019], Lemma 4, Lemma 5, Proposition 3\]は，$\{ \mu_k \}_{k \ge k(\hat{\ell})}$が$[\underline{\mu}, \bar{\mu}]$内の非減少列であることを仮定している．null step後に$\mu_k$を減らすことが禁止されているのは，まさにこの仮定を成り立たせるためである．
:::

-->

<!--


### $h$が微分可能凸関数の最大の場合：d-stationary点を求める版

第2成分$h$が，有限個の**微分可能**な凸関数の最大

$$
    h(x) = \max_{1 \le i \le N} \psi_i(x)
$$

で書けている場合には，臨界点よりも強いd-stationary性を保証できる．このとき$h$の劣微分は活性な勾配の凸包

$$
    \partial h(x) = \mathrm{conv} \left( \{ \nabla \psi_i(x) \}_{i \in A(x)} \right) , \quad A(x) := \{ 1 \le i \le N : \psi_i(x) = h(x) \}
$$

である．したがって，すべての$i \in A(\bar{x})$について

$$
    \bar{x} \in \argmin_{x \in X} \left( g(x) - \langle \nabla \psi_i(\bar{x}) , x \rangle \right)
$$

が成り立てば，$\bar{x}$はd-stationary点になる．DCAの1ステップが，活性なすべての$i$について同時に$\bar{x}$を返す，という条件である．そこで，1つの劣勾配$\nabla h(x^{k(\ell)})$を使う代わりに，**活性なすべての$\nabla \psi_i(x^{k(\ell)})$について部分問題を解く**．技術的な理由から，活性集合$A(x)$そのものではなく，$\epsilon > 0$による緩和

$$
    A_\epsilon(x) := \{ 1 \le i \le N : \psi_i(x) \ge h(x) - \epsilon \}
$$

を使う（$\epsilon = 0$では，d-stationaryでない臨界点に留まりうる）．$i \in A_\epsilon(x^{k(\ell)})$ごとの部分問題は

$$
    \min_{x \in X, \, r \in \mathbb{R}} \ r - \langle \nabla \psi_i(x^{k(\ell)}) , x - x^{k(\ell)} \rangle + \frac{\mu_k}{2} \| x - x^{k(\ell)} \|_2^2
    \quad \text{s.t.} \quad g(x^j) + \langle \nabla g(x^j) , x - x^j \rangle \le r , \ \ j \in B^k
    \tag{\mathrm{P}_i}
$$

であり，$(P)$の$\nabla h(x^{k(\ell)})$を$\nabla \psi_i(x^{k(\ell)})$で置き換えただけのものである．

#### アルゴリズム2

**入力**：$x^0 \in X$，$\kappa \in (0,1)$，$0 < \underline{\mu} \le \mu_0 \le \bar{\mu} < \infty$，許容誤差$\delta_{\mathrm{Tol}} \ge 0$，緩和パラメータ$\epsilon > 0$

1. $g(x^0), \, \nabla g(x^0), \, h(x^0)$と，すべての$i \in A_\epsilon(x^0)$に対する$\nabla \psi_i(x^0)$を計算
1. $B^0 \gets \{ 0 \}$，$k(0) \gets \ell \gets 0$
1. for $k = 0, 1, 2, \ldots$
1. 　　for all $i \in A_\epsilon(x^{k(\ell)})$
1. 　　　　$(\mathrm{P}_i)$を解いて，一意解$(r_i , y_i)$と乗数$\alpha^i \in \mathbb{R}_+^{|B^k|}$を得る
1. 　　endfor
1. 　　$x^{k+1} \gets y_{i^*}$，ただし $i^* \in \argmin_{i \in A_\epsilon(x^{k(\ell)})} \left[ r_i - \langle \nabla \psi_i(x^{k(\ell)}) , y_i - x^{k(\ell)} \rangle + \dfrac{\mu_k}{2} \| y_i - x^{k(\ell)} \|_2^2 \right]$
1. 　　**［停止判定］** if $\| x^{k+1} - x^{k(\ell)} \|_2 \le \delta_{\mathrm{Tol}}$ then return $x^{k(\ell)}$
1. 　　**［オラクル］** $g(x^{k+1})$，$\nabla g(x^{k+1})$，$h(x^{k+1})$を計算
1. 　　if $f(x^{k+1}) \le f(x^{k(\ell)}) - \dfrac{\kappa \underline{\mu}}{2} \| x^{k+1} - x^{k(\ell)} \|_2^2$
1. 　　　　**［serious step］** すべての$i \in A_\epsilon(x^{k+1})$に対し$\nabla \psi_i(x^{k+1})$を計算
1. 　　　　$x^{k(\ell+1)} \gets x^{k+1}$，$k(\ell+1) \gets k+1$，$\ell \gets \ell + 1$
1. 　　　　$B^{k+1} \supset \{ k+1 \}$と$\mu_{k+1} \in [ \underline{\mu} , \bar{\mu} ]$を選ぶ
1. 　　else
1. 　　　　**［null step］** $B^{k+1} \supset \{ -k , \, k(\ell) , \, k+1 \}$と$\mu_{k+1} \in [ \mu_k , \bar{\mu} ]$を選ぶ
1. 　　endif
1. endfor

#### アルゴリズム1との違い

- **1反復あたり$|A_\epsilon(x^{k(\ell)})|$個の部分問題を解く**（高々$N$個）．$N = 1$なら，降下判定を除いてアルゴリズム1と一致する．
- **降下判定が$f$そのものになっている**（10行目）．$\epsilon > 0$のため$\nabla \psi_{i^*}(x^{k(\ell)})$が$\partial h(x^{k(\ell)})$に属さないことがありえて，そのときアルゴリズム1の$g$だけの判定からは$f$の減少が導けないからである．代償として，**$h$の値を毎反復評価する必要がある**．
- **保証が強い**：$\delta_{\mathrm{Tol}} = 0$かつ準位集合が有界なら，安定中心の列の任意の集積点は**d-stationary点**である．
- **必要なオラクルも強い**：$\partial h$の元を1つ返す$\nabla h$では足りず，活性な$\psi_i$の勾配をすべて列挙できる必要がある．すなわち，記事の冒頭で置いた問題設定より強い仮定を課している．
- serious step後の$\mu_{k+1}$に，アルゴリズム1のような上限$\mu_k$（＝減少のみ）は課されていない．

なお，null stepでバンドルに加える集約線形化$\bar{g}^{-k}$は，複数ある部分問題のうち採用された$i^*$の乗数$\alpha^{i^*}$から作る（$p^{k+1} = \sum_{j \in B^k} \alpha^{i^*}_j \nabla g(x^j)$）．論文はこの点を明示していないが，$x^{k+1} = y_{i^*}$が実際に解いた問題の最適性条件から来るのはこの組み合わせだけである．

-->

## 数値実験

\[[dO2019][dO2019], Sect. 6\]の要約である．数値はすべて同論文の表・図から引いたもので，各主張の出所を都度示す．以下，アルゴリズム1の実装を**PBM1**と呼ぶ．

:::message
この節では，適宜\[[dO2019][dO2019]\]の図や表を参照してほしい．
:::

### 設定

比較対象は次の8つである（\[[dO2019][dO2019], Sect. 6\]）．

| ソルバ | 中身 |
| --- | --- |
| PBM1 | アルゴリズム1（$g$を切除平面モデル$\check{g}^k$で，$h$を線形化1本で近似する．部分問題$(P)$は**凸**） |
| PBM3 | アルゴリズム1の変種で，$h$を3本の線形化で近似する．部分問題は**非凸**になる |
| DCA-CPM / DCA-LBM | 素のDCAの部分問題$\argmin_{x \in X}(g(x) - \langle \nabla h(x^k), x\rangle)$を，それぞれKelleyの切除平面法／level bundle法で解く実装 |
| PLM | 近接線形化法（DCAに近接項を足したもの，部分問題は切除平面近似で解く） |
| HANSO | DC構造を使わない一般の非平滑非凸ソルバ（BFGS＋gradient sampling） |
| PBDC / DCPCA | 他文献のバンドル法．PBM3と同じく$h$を複数の線形化で近似する |

実装はMATLAB 2017a＋Gurobi 7.5.1，停止判定は$\delta_{\mathrm{Tol}} = 10^{-4}$，時間制限は1問あたり3600秒である．バンドルの上限本数は$\max\{ 100, \min\{ d+5, 1000 \} \}$で，serious stepの後は活性なカットのみを残す（\[[dO2019][dO2019], Sect. 6\]）．

:::message alert
論文の実装は，本記事で説明した降下判定$(D)$（\[[dO2019][dO2019], 式 (13)\]）ではなく，$f$そのものを使う判定

$$ f(x^{k+1}) \le f(x^{k(\ell)}) - \frac{\kappa \underline{\mu}}{2} \| x^{k+1} - x^{k(\ell)} \|_2^2 $$

（\[[dO2019][dO2019], 式 (12)\]，$\kappa = 0.1$，$\underline{\mu} = 10^{-5}$）を採用している．したがって数値実験では「$h$の値を一度も評価しない」というアルゴリズム1の利点は使われておらず，$g$と$h$の評価回数は同じになっている．
:::

### 問題群

| 問題群 | 内容 | 次元$d$ |
| --- | --- | --- |
| 無制約DC<br>\[[dO2019][dO2019], Sect. 6.1, Table 1\] | 先行研究で標準的に使われる10問 | 2〜50,000 |
| 線形制約DC<br>\[[dO2019][dO2019], Sect. 6.2, Table 2, Table 3, Table 4\] | 確率制約問題$\mathbb{P}[c(x,\xi) \le 0] \ge p$のDC近似．期待値はMonte-Carlo（10,000シナリオ）で近似．エネルギー計画問題PlanToyとノルム最適化問題，各18インスタンス | 2〜200 |
| 2次制約DC<br>\[[dO2019][dO2019], Sect. 6.3, Table 5\] | $\frac{1}{2}\sum_i \alpha_i x_i^2 \le \frac{1}{2}K^2$を制約に持つDC計画．$(P)$がQCQPになる | 2〜50,000 |

3番目は$X$が多面体でない例で，$(P)$がQPでなくQCQPになる．本手法が$X$の構造に応じた専用ソルバを使える設計になっていることの実証になっている．

### 結果

**無制約（\[[dO2019][dO2019], Sect. 6.1\]）．** performance profile（\[[dO2019][dO2019], Fig. 1\]）で見ると，

- $\nabla g$の評価回数と頑健性では，4つのバンドル法（PBM1, PBM3, PBDC, DCPCA）が他を上回る
- $\nabla h$の評価回数では，逆にDCA-CPM／DCA-LBMが上回る（DCAは1反復に$\nabla h$を1回しか使わないため）
- CPU時間と頑健性のバランスではPBM1が最も良い

また\[[dO2019][dO2019], Table 1\]の問題7・8・9では，DCA系とPLMが既知の最良値に届かなかったのに対し，PBM1は届いた．

**線形制約（\[[dO2019][dO2019], Sect. 6.2\]）．** PlanToyの18インスタンスすべてで（\[[dO2019][dO2019], Table 2\]），PLMは大域解に到達できないだけでなく，**実行可能ですらない点**で停止した（罰金付き定式化の臨界点に落ちるため）．PBM1・PBM3は成功し，得られた点は確率制約をほぼ満たしていた．ノルム最適化問題（\[[dO2019][dO2019], Table 3, Table 4\]）では，DCA-CPMは$d \ge 50$で臨界点の計算に失敗，DCA-LBMとPLMは臨界点は出すが大域解には届かず，HANSOは関数評価が多すぎて$d \ge 100$では1時間以内に終わらなかった．PBM1は全インスタンスで既知の最適値に近い値を得ている．

**2次制約（\[[dO2019][dO2019], Sect. 6.3, Table 5\]）．** PBM1が最速である．理由は1反復あたりに解く部分問題の個数で，PBM1は$\#\nabla g$個のQCQP，PBM3はその3倍，DCA-CPMは$\#\nabla h$個のLP，PLMは$\#\nabla h$個のQCQPを解く．PLMは$d \ge 100$で1時間以内に解けていない．

### 読み取れること

- **$h$の線形化は1本で十分**：PBM1（部分問題は凸QP 1個）は，$h$を複数の線形化で近似するPBM3・PBDC・DCPCAに対して反復回数でも計算時間でも劣らず，解の質も落ちない．「$h$は安定中心での線形化1本でよい」という設計判断が実験で裏付けられている．
- **DC構造を使う利点**：DC構造を使わないHANSOは，関数・劣勾配の評価回数が桁違いに多い（\[[dO2019][dO2019], Table 1\]の問題4，$d = 100$で$45{,}619$回，対するPBM1は$101$回）．
- **大域最適性の保証はない**：\[[dO2019][dO2019], Table 1\]の$\bar{f}$は既知の最良値で，$*$印は到達できなかったことを示す．PBM1にも$*$のつく例はある．逆に，一般の非凸ソルバHANSOだけが大域解に到達した問題も1つある（\[[dO2019][dO2019], Sect. 6.1\]）．定理 3. が保証するのは臨界性までであり，実験結果もその範囲を超えない．

なお\[[dO2019][dO2019], Sect. 6.4\]は，$h$が微分可能凸関数の最大で書ける問題に対する，同論文のAlgorithm 2（PBM-d．臨界点より強いd-stationary点を求める変種で，本記事では扱わない）の実験である．臨界点だがd-stationaryでない点$\tilde{x} = (0,0)$を持つ2次元問題で，PBM1は8つの初期点のうちいくつかで$\tilde{x}$に停止するのに対し，PBM-dは$\tilde{x}$の近くで2つの部分問題を解くため脱出できる（\[[dO2019][dO2019], Fig. 2\]）．

## 関連研究

### DCAからの系譜

本手法は，DCAに2段階の修正を加えたものと見ると位置づけがはっきりする．

| 手法 | 部分問題 |
| --- | --- |
| DCA \[TL1997\]（サーベイは\[[LP2018][LP2018]\]） | $\argmin_{x \in X} \left( g(x) - \langle \nabla h(x^k) , x \rangle \right)$ |
| 近接線形化法（PLM）\[[SOS2016][SOS2016]\] | 上に近接項$\mu_k D(x, x^k)$を追加．\[[CNOSS2018][CNOSS2018]\]が一般のBregman関数$D$に拡張 |
| \[[dO2019][dO2019]\] アルゴリズム1 | さらに$g$を切除平面モデル$\check{g}^k$に置き換え，中心を安定中心$x^{k(\ell)}$に固定 |

PLMまでは，毎反復で非平滑凸計画を**厳密に**解く必要が残る．\[[SOS2016][SOS2016]\]は非厳密版を扱っているが，部分問題を漸近的に最適まで解くことに加え，計算される$g$の劣勾配が直前の$h$の劣勾配に対して近接性条件を満たすことを要求する．\[[dO2019][dO2019]\]のアルゴリズム1は，オラクルが返す$\nabla g(x^{k+1})$に**何の条件も課さない**点が違う．切除平面モデルを使うことで，部分問題を近似的に解くという操作が「$g$をモデルで置き換えて厳密に解く」という形に整理され，近似の度合いを制御する必要がなくなっている．

### $h$を複数の線形化で近似する手法との比較

$g$だけでなく$h$の側も複数の線形化で近似する，すなわち$\check{h}^k := \max_{j \in B_h^k} \bar{h}^j$（$\lvert B_h^k \rvert \ge 2$）を使って$f$を$\check{g}^k - \check{h}^k$で近似するバンドル法が先行している．数値実験に現れたPBDC \[[JBKM2017][JBKM2017]\]とDCPCA \[[GGMB2018][GGMB2018]\]がそれである．

$-\check{h}^k = \min_j (-\bar{h}^j)$は凹の折れ線なので，$\check{g}^k - \check{h}^k$は真に非凸になり，その大域最小化には枝ごとの列挙が要る．一方アルゴリズム1では$h$の近似$\bar{h}^{k(\ell)}$がアフィンなので，$\check{g}^k - \bar{h}^{k(\ell)}$は凸のままである（凸関数からアフィン関数を引いても凸性は保たれる）．この違いが以下の差をすべて生んでいる．

| | アルゴリズム1 | PBDC \[[JBKM2017][JBKM2017]\] | DCPCA \[[GGMB2018][GGMB2018]\] |
| --- | --- | --- | --- |
| 実行可能集合 | 凸集合$X$ | $X = \mathbb{R}^d$のみ | $X = \mathbb{R}^d$のみ |
| $h$の近似 | 線形化1本（アフィン） | 線形化複数本（折れ線） | 線形化複数本（折れ線） |
| 部分問題の凸性 | 凸 | 非凸 | 非凸 |
| 1反復の部分問題 | 狭義凸計画1個 | QP $\lvert B_h^k \rvert$個（最低2個） | QP 2個＋直線探索 |
| 部分問題の型 | $X$上の狭義凸計画（QP，QCQP，錐計画など） | QP | QP |
| 直線探索 | 不要 | 不要 | 必要 |
| $h$の値の評価 | 不要 | 必要 | 必要 |

アルゴリズム1の利点は，$f$への近似の精度ではなく適用範囲と1反復あたりのコストにある．

1. **制約付き問題を扱える．** PBDC・DCPCAは無制約DC計画専用である．本記事の設定（$X$は有界凸集合）ではそもそも適用できないので，これが決定的な差になる．
2. **1反復に解く部分問題が1個．** 上に述べたとおり$\check{g}^k - \check{h}^k$は非凸なので，PBDCは$\lvert B_h^k \rvert$個のQPを解いて大域最小化する．近似を良くするほど個数が増える．$h$を線形化1本に固定する割り切りが，この列挙を消している．
3. **部分問題がQPでなくてよい．** 近接項のBregman関数$\omega$を自由に選べるため，$(P)$は「$X$上の狭義凸計画」であればよい．$X$が球・単体・スペクトラヘドロンなどのときは専用ソルバを使える（\[[dO2019][dO2019], Sect. 1\]）．実際\[[dO2019][dO2019], Sect. 6.3\]では$(P)$がQCQPになる問題を解いている．
4. **直線探索もLipschitz定数の推定も不要．**
5. **$h$の値を評価しない．** $\check{h}^k$を作るには$h$の値と複数の劣勾配が要る．

代償は近似の粗さである．$\check{g}^k - \bar{h}^{k(\ell)}$は$\check{g}^k - \check{h}^k$より$f$の近似としては悪いはずで，反復回数が増えても不思議はない．
すなわち$h$の近似を粗くしても反復回数・計算時間・解の質は悪化しなかった，という実験結果をもって，この単純化を正当化している．ただし全面的な優位ではなく，DCPCAは$\nabla g$の評価回数ではPBM1より少ない（代わりに$\nabla h$が多い）（\[[dO2019][dO2019], Sect. 6.1\]）．

### 停留性の強さ

理論保証の面では，PBDCも臨界点までで定理 3. と同格である．より強い停留性を狙う研究として次がある．

- \[[JBKMT2018][JBKMT2018]\]：無制約DC計画に対しClarke停留点を保証する．ただし任意の点で$\partial g, \partial h$が多面体であることを仮定する（$g, h$がともに有限個の微分可能関数の最大である場合などが該当）．
- \[[PRA2017][PRA2017]\]：d-stationary点を求めるアルゴリズムを提案した．\[[dO2019][dO2019]\]のAlgorithm 2はこれに触発されたものである．

停留性の概念には 局所最小解 $\Rightarrow$ d-stationary $\Rightarrow$ Clarke stationary $\Rightarrow$ 臨界点 という強さの順序がある（\[[dO2019][dO2019]\] 2節）．

### DC構造を使わない手法

数値実験の比較対象HANSO \[[LO2013][LO2013]\]は，BFGSとgradient samplingによる一般の非平滑非凸最適化ソルバで，DC分解を一切使わない．適用範囲は広いが，関数・劣勾配の評価回数が桁違いに多くなる（\[[dO2019][dO2019], Table 1\]）．目的関数の評価が高価な問題では実用にならない一方，DC構造を使う手法が臨界点に留まる場面で大域解に到達できた例もある（\[[dO2019][dO2019], Sect. 6.1\]）．

# 参考文献

\[[dO2019][dO2019]\] W. de Oliveira, Proximal bundle methods for nonsmooth DC programming, Journal of Global Optimization, 75(2):523–563, 2019年．

\[[CNOSS2018][CNOSS2018]\] J. X. Cruz Neto, P. R. Oliveira, A. Soubeyran, J. C. O. Souza, A generalized proximal linearized algorithm for DC functions with application to the optimal size of the firm problem, Annals of Operations Research, 2018年．

\[[GGMB2018][GGMB2018]\] M. Gaudioso, G. Giallombardo, G. Miglionico, A. M. Bagirov, Minimizing nonsmooth DC functions via successive DC piecewise-affine approximations, Journal of Global Optimization, 71(1):37–55, 2018年．（DCPCA）

\[[JBKM2017][JBKM2017]\] K. Joki, A. M. Bagirov, N. Karmitsa, M. M. Mäkelä, A proximal bundle method for nonsmooth DC optimization utilizing nonconvex cutting planes, Journal of Global Optimization, 68(3):501–535, 2017年．（PBDC）

\[[JBKMT2018][JBKMT2018]\] K. Joki, A. Bagirov, N. Karmitsa, M. M. Mäkelä, S. Taheri, Double bundle method for finding Clarke stationary points in nonsmooth DC programming, SIAM Journal on Optimization, 28(2):1892–1919, 2018年．

\[[LP2018][LP2018]\] H. A. Le Thi, T. Pham Dinh, DC programming and DCA: thirty years of developments, Mathematical Programming, 169(1):5–68, 2018年．（DCAの30年をまとめたサーベイ）

\[[LO2013][LO2013]\] A. S. Lewis, M. L. Overton, Nonsmooth optimization via quasi-Newton methods, Mathematical Programming, 141(1):135–163, 2013年．（HANSO）

\[[PRA2017][PRA2017]\] J.-S. Pang, M. Razaviyayn, A. Alvarado, Computing B-stationary points of nonsmooth DC programs, Mathematics of Operations Research, 42(1):95–118, 2017年．

\[[SOS2016][SOS2016]\] J. C. O. Souza, P. R. Oliveira, A. Soubeyran, Global convergence of a proximal linearized algorithm for difference of convex functions, Optimization Letters, 10(7):1529–1539, 2016年．

\[TL1997\] P. D. Tao, H. A. Le Thi, Convex analysis approach to DC programming: theory, algorithms and applications, Acta Mathematica Vietnamica, 22(1):289–355, 1997年．

\[[Tri2023][Tri2023]\] Trigger-FK, 初学者が学んだDC最適化（概要）, 2023年．

[dO2019]:https://doi.org/10.1007/s10898-019-00755-4
[CNOSS2018]:https://doi.org/10.1007/s10479-018-3104-8
[GGMB2018]:https://doi.org/10.1007/s10898-017-0568-z
[JBKM2017]:https://doi.org/10.1007/s10898-016-0488-3
[JBKMT2018]:https://doi.org/10.1137/16M1115733
[LO2013]:https://doi.org/10.1007/s10107-012-0514-2
[LP2018]:https://doi.org/10.1007/s10107-018-1235-y
[PRA2017]:https://doi.org/10.1287/moor.2016.0795
[SOS2016]:https://doi.org/10.1007/s11590-015-0969-1
[Tri2023]:https://qiita.com/Trigger-FK/items/84148dcf2c92e9647485



<!---
凸関数f,gのオラクルと，その劣勾配\nabla f, \nabla g のオラクルは計算できるものとします．
このとき，これらのオラクルと，領域Xの性質を使って\min_{x \in X} ( f(x)-g(x))を計算するようなDCAアルゴリズムとして，代表的な手法および，現実的に速い方法を列挙してください


まず重要なのは、\(F(x)=f(x)-g(x)\) は一般には非凸なので、**通常のDCAで大域的最適解を計算できるとは限らない**という点です。標準DCAが保証するのは、主としてDC臨界点、すなわち
\[
\partial f(\bar x)\cap \partial g(\bar x)\neq\varnothing
\]
を満たす点への収束です。DC分解や制約集合の構造によって、実用上の速度と大域性は大きく変わります。DC最適化の標準的な定式化とDCAの基本的性質については、例えば田中・奥野の概説が参考になります。 [jstage.jst.go](https://www.jstage.jst.go.jp/article/bjsiam/29/3/29_14/_pdf/-char/ja)

以下では、\(X\subseteq\mathbb{R}^n\) は閉凸集合とし、制約を
\[
\iota_X(x)=
\begin{cases}
0,&x\in X,\\
+\infty,&x\notin X
\end{cases}
\]
に吸収して
\[
\min_x\; f(x)+\iota_X(x)-g(x)
\]
と書きます。

## 1. 標準DCA

現在点 \(x_k\in X\) で
\[
y_k\in \partial g(x_k)
\]
を計算し、次の凸問題を解きます。
\[
x_{k+1}\in\arg\min_{x\in X}
\left\{f(x)-\langle y_k,x-x_k\rangle\right\}.
\tag{DCA}
\]

定数項を除けば、
\[
x_{k+1}\in\arg\min_{x\in X}
\{f(x)-\langle y_k,x\rangle\}.
\]

これは、凸関数 \(g\) を接平面
\[
g(x_k)+\langle y_k,x-x_k\rangle
\]
で下から線形近似し、得られた凸上界を最小化するMM法と解釈できます。目的関数値は、子問題を厳密に解けば
\[
F(x_{k+1})\leq F(x_k)
\]
となります。 [jstage.jst.go](https://www.jstage.jst.go.jp/article/bjsiam/29/3/29_14/_pdf/-char/ja)

### 擬似コード

```text
x0 ∈ X を選ぶ
for k = 0,1,2,...
    yk ∈ ∂g(xk) を計算
    xk+1 ≈ argmin_{x∈X} f(x) - <yk, x>
    if ||xk+1 - xk|| が十分小さい:
        return xk+1
end
```

停止条件としては、例えば
\[
\frac{\|x_{k+1}-x_k\|}{\max\{1,\|x_k\|\}}\leq \varepsilon_x
\]
に加えて、
\[
|F(x_{k+1})-F(x_k)|\leq \varepsilon_F
\]
を用いるのが自然です。

ただし、標準DCAでは各反復で凸子問題を十分正確に解く必要があります。大規模問題では、この子問題の解法がボトルネックになることが多く、子問題を解析的または近接写像で処理できるDC分解を選ぶことが重要です。 [jstage.jst.go](https://www.jstage.jst.go.jp/article/bjsiam/29/3/29_14/_pdf/-char/ja)

## 2. Fenchel双対型DCA

標準DCAは共役関数を使うと
\[
y_k\in\partial g(x_k),\qquad
x_{k+1}\in\partial f^*(y_k)
\]
と書けます。後半は
\[
x_{k+1}\in\arg\min_{x\in X}\{f(x)-\langle y_k,x\rangle\}
\]
と同値です。

したがって、以下の2つのオラクルが効率的なら非常に自然です。

- \(g\) の劣勾配オラクル：\(y_k\in\partial g(x_k)\)
- \(f+\iota_X\) の線形摂動最小化オラクル：
  \[
  \arg\min_{x\in X}\{f(x)-\langle y_k,x\rangle\}
  \]

特に \(f\) が線形計画、二次計画、エントロピー型関数などで、線形項を加えた最小化が高速に解ける場合に有効です。凸共役を陽に計算する必要はなく、子問題を直接解けばよいという点が実装上重要です。 [jstage.jst.go](https://www.jstage.jst.go.jp/article/bjsiam/29/3/29_14/_pdf/-char/ja)

## 3. 近接DCA

標準DCAの子問題が重い場合、二次正則化を加えます。
\[
x_{k+1}\approx
\arg\min_{x\in X}
\left\{
f(x)-\langle y_k,x-x_k\rangle
+\frac{\rho_k}{2}\|x-x_k\|^2
\right\},
\qquad
y_k\in\partial g(x_k).
\tag{pDCA}
\]

\(\rho_k>0\) により子問題は強凸になり、解が一意になりやすく、反復解法も安定します。

例えば \(f\) が微分可能で \(\nabla f\) が \(L\)-Lipschitz なら、
\[
x_{k+1}
=
\operatorname{prox}_{(f+\iota_X)/\rho_k}
\left(x_k+\frac{y_k}{\rho_k}\right)
\]
という形になります。さらに \(f\) の勾配を線形化する場合には、
\[
x_{k+1}
\in\arg\min_{x\in X}
\left\{
\langle \nabla f(x_k)-y_k,x-x_k\rangle
+\frac{L}{2}\|x-x_k\|^2
\right\},
\]
すなわち
\[
x_{k+1}
=
\operatorname{prox}_{\iota_X/L}
\left(x_k-\frac{1}{L}(\nabla f(x_k)-y_k)\right)
\]
となります。

\(X=\mathbb{R}^n\) なら単に
\[
x_{k+1}=x_k-\frac{1}{L}(\nabla f(x_k)-y_k)
\]
です。\(X\) が単純集合なら、射影
\[
x_{k+1}
=
\Pi_X\left(x_k-\frac{1}{L}(\nabla f(x_k)-y_k)\right)
\]
で計算できます。

## 4. 外挿付きDCA・BDCA

反復が遅い場合には、外挿を加えます。典型的には
\[
z_k=x_k+\beta_k(x_k-x_{k-1}),
\qquad 0\leq \beta_k<1,
\]
で外挿点を作り、
\[
y_k\in\partial g(z_k)
\]
または \(y_k\in\partial g(x_k)\) を用いて凸子問題を解きます。

滑らかな \(f\) と近接可能な凸部分 \(h\) を含む問題
\[
\min_x f(x)+h(x)-g(x)
\]
に対するpDCAeの典型形は
\[
\begin{aligned}
z_k&=x_k+\beta_k(x_k-x_{k-1}),\\
\xi_k&\in\partial g(x_k),\\
x_{k+1}
&=
\arg\min_x
\left\{
\langle\nabla f(z_k)-\xi_k,x\rangle
+\frac{L}{2}\|x-z_k\|^2+h(x)
\right\}.
\end{aligned}
\]

外挿の計算コストは標準pDCAとほぼ同じですが、実際の数値性能が改善することがあります。ただし、非凸問題では無条件のNesterov加速は不安定になり得るため、目的関数の増加時に \(\beta_k=0\) に戻すリスタートを併用するのが実用的です。 [orsj](https://orsj.org/wp-content/corsj/or69-3/or69_3_136.pdf)

また、より単純な方法として、標準DCAで得た方向
\[
d_k=y_k^{\mathrm{DCA}}-x_k
\]
に対して
\[
x_{k+1}=y_k^{\mathrm{DCA}}+\lambda_k d_k
\]
と進み、目的関数減少条件を満たすまで \(\lambda_k\) をバックトラックするBDCA型の方法もあります。

## 5. 非厳密DCA

子問題を厳密に解かず、
\[
x_{k+1}\approx\arg\min_{x\in X}
\{f(x)-\langle y_k,x\rangle\}
\]
としても、誤差を適切に減少させれば収束性を維持できます。

例えば、子問題の最適性残差に対して
\[
\operatorname{dist}\left(
0,\partial(f+\iota_X)(x_{k+1})-y_k
\right)
\leq \eta_k
\]
を課し、
\[
\eta_k\to0
\quad\text{または}\quad
\sum_{k=0}^\infty \eta_k<+\infty
\]
とします。

大規模問題では、各外側反復で子問題を完全に解くよりも、

- 内側のAPG/FISTA
- 投影勾配法
- 座標降下法
- ADMM
- 半滑らかNewton法
- L-BFGS型の近似Newton法

を数回だけ実行する方が、全体として速いことが多いです。近年の解説でも、非厳密Newton型の二次近似子問題と直線探索を組み合わせる方法が扱われています。 [orsj](https://orsj.org/wp-content/corsj/or69-3/or69_3_136.pdf)

## 6. DC分解を変更して子問題を簡単にする

これは実用上、最も重要な工夫の一つです。DC分解は一意ではなく、任意の凸関数 \(\phi\) に対して
\[
f-g=(f+\phi)-(g+\phi)
\]
も同じ目的関数のDC分解です。

例えば \(f\) が微分可能で \(\nabla f\) が \(L\)-Lipschitz なら、
\[
f(x)
=
\frac{L}{2}\|x\|^2
-
\left(
\frac{L}{2}\|x\|^2-f(x)
\right)
\]
という分解が使えます。後者も凸になります。

この分解では、\(\partial g(x_k)\) に相当する項が
\[
Lx_k-\nabla f(x_k)
\]
となるので、DCA更新は
\[
x_{k+1}
\in
\arg\min_{x\in X}
\left\{
f_{\mathrm{simple}}(x)
+\frac{L}{2}\|x-x_k\|^2
+\langle \nabla f(x_k)-\xi_k,x-x_k\rangle
\right\}
\]
のような近接勾配型になります。

代表例は次の通りです。

| \(X\)・正則化の構造 | 典型的な子問題解法 |
|---|---|
| \(X=\mathbb{R}^n\)、\(\ell_1\)項 | soft-thresholding |
| 箱制約 \(l\leq x\leq u\) | clipping付きsoft-thresholding |
| 単体 \(X=\{x\geq0,\mathbf{1}^\top x=1\}\) | simplex projection |
| \(\ell_2\)球 | 球への射影 |
| affine制約 | 線形方程式または射影 |
| separable \(f\) | 座標別更新 |
| sparse/group-sparse構造 | block soft-thresholding |
| 近接写像が既知 | prox評価1回 |

例えば \(X=\mathbb{R}^n\) で
\[
f(x)=\ell(x)+\lambda\|x\|_1,\qquad
g(x)\ \text{が凸}
\]
なら、近接DCAは
\[
x_{k+1}
=
\operatorname{soft}_{\lambda/L}
\left(
x_k-\frac{1}{L}
\left(\nabla\ell(x_k)-y_k\right)
\right)
\]
という形にできます。このように、**一般の凸ソルバを毎回呼ぶのではなく、射影・近接写像だけで更新する設計**が、現実には最も高速になりやすいです。 [jstage.jst.go](https://www.jstage.jst.go.jp/article/bjsiam/29/3/29_14/_pdf/-char/ja)

## 7. 大規模有限和構造なら確率的DCA

もし
\[
f(x)=\frac1N\sum_{i=1}^N f_i(x),
\qquad
g(x)=\frac1N\sum_{i=1}^N g_i(x)
\]
のような有限和構造があるなら、毎回すべての劣勾配を計算する必要はありません。

代表的には、

- stochastic DCA
- DCA-SVRG
- DCA-SAGA
- DCA-SAG
- minibatch DCA

があります。

例えば \(g\) の劣勾配をSVRG型に近似して
\[
y_k
\approx
\frac1b\sum_{i\in I_k}
\bigl(
\partial g_i(x_k)-\partial g_i(\widetilde x)
\bigr)
+\partial g(\widetilde x)
\]
とし、
\[
x_{k+1}
\approx
\arg\min_{x\in X}
\{f(x)-\langle y_k,x\rangle\}
\]
を解きます。

DCA-SVRGやDCA-SAGAでは、適切な強凸性・Lipschitz条件の下でDC臨界点への収束が解析されています。大規模有限和問題では、フルバッチDCAより実用的になる可能性があります。 [jmlr](https://jmlr.org/papers/volume23/21-1146/21-1146.pdf)

## 8. 多面体的構造がある場合

\(g\) が有限個のアフィン関数の最大
\[
g(x)=\max_{1\leq i\leq m}\{a_i^\top x-b_i\}
\]
で表される場合、
\[
\min_x f(x)-g(x)
=
\min_x\min_i
\{f(x)-a_i^\top x+b_i\}
\]
なので、
\[
\min_x f(x)-g(x)
=
\min_i
\left[
\min_{x\in X}\{f(x)-a_i^\top x+b_i\}
\right]
\]
と分解できます。

したがって、有限個の凸問題を解けば大域的最適解が得られます。これは標準DCAの局所的な性質とは異なり、特定の多面体構造に基づく大域解法です。ただし、\(m\) が大きい場合には全列挙が高価なので、列生成、アクティブ集合、branch-and-bound、cutting-planeを組み合わせます。 [jstage.jst.go](https://www.jstage.jst.go.jp/article/bjsiam/29/3/29_14/_pdf/-char/ja)

## 9. 大域的最適解が必要な場合

一般の凸関数オラクルと劣勾配オラクルだけから、任意のDC問題の大域的最適解を効率的に得ることは期待できません。DCAの出力は基本的に、初期点とDC分解に依存する局所解・臨界点です。

大域解を狙う代表的な方法は次の通りです。

1. **多点初期化DCA**  
   複数の初期点からDCAを実行し、最良値を採用する。

2. **DCA + 局所探索**  
   DCAで候補点を得た後、bundle法、trust-region法、SQP型法などで改善する。

3. **分枝限定法**  
   \(g\) の凸下界、\(f-g\) の領域ごとの下界を用いて枝刈りする。

4. **凸緩和・cutting-plane**  
   \(g\) を有限個の接平面で近似し、下界問題を解く。

5. **多面体DC専用法**  
   有限個の凸子問題、列生成、アクティブなアフィン片の探索を利用する。

6. **凹最小化型の特殊アルゴリズム**  
   \(f\) が線形または単純な凸関数で、\(-g\) の凹性が支配的な場合に利用する。

## 10. 実務上の推奨順位

与えられたオラクルが \(f,g,\partial f,\partial g\) と、\(X\) に関する射影・線形最小化・近接写像のいずれかであるとすると、次の順で試すのが現実的です。

### ケースA：近接写像または射影が高速

\[
\boxed{\text{近接DCAまたはpDCAe}}
\]

- \(\nabla f\) または \(\partial f\) を使う。
- 二次正則化を入れる。
- 外挿とリスタートを使う。
- 子問題を1回のprox・射影で処理する。

最も実装しやすく、大規模問題に向いています。

### ケースB：線形摂動付き凸最小化が高速

\[
\boxed{\text{標準DCA}}
\]

\[
x_{k+1}\in\arg\min_{x\in X}
\{f(x)-\langle y_k,x\rangle\}.
\]

LP、QP、エントロピー最小化、ネットワークフローなど、子問題を高速に解ける場合はこちらが有力です。

### ケースC：有限和構造が大きい

\[
\boxed{\text{DCA-SVRG / DCA-SAGA}}
\]

フル勾配・フル劣勾配の計算が支配的なら、variance reduction付きの確率的DCAを使います。

### ケースD：子問題が高価だが二次近似が利用可能

\[
\boxed{\text{inexact proximal DCA + line search}}
\]

子問題を完全には解かず、内側反復を少数回で打ち切ります。各外側反復のコストを抑えられます。

### ケースE：大域解が必要

\[
\boxed{\text{多点DCAを初期候補にし、branch-and-boundまたはcutting-planeを併用}}
\]

通常のDCA単独を大域解法とみなしてはいけません。

要するに、最も標準的な更新は
\[
y_k\in\partial g(x_k),\qquad
x_{k+1}\in\arg\min_{x\in X}
\{f(x)-\langle y_k,x\rangle\},
\]
最も実用的な高速化は
\[
x_{k+1}
=
\operatorname{prox}_{\text{simple}/L}
\left(
x_k-\frac{1}{L}\bigl(\nabla f(x_k)-y_k\bigr)
\right)
\]
型に落とし込むことです。したがって、\(X\) の射影、\(f\) の近接写像、あるいは線形摂動付きの凸最小化のどれが高速に計算できるかを基準に、標準DCA・pDCAe・非厳密DCAを選ぶのがよいでしょう。
--->

<!---
f,g はLipschitzまでしか仮定できません.　C^2級は仮定できません．

はい。その場合、前回の回答にあった

\[
x_{k+1}
=
\Pi_X\left(x_k-\frac{1}{L}\bigl(\nabla f(x_k)-y_k\bigr)\right)
\]

のような更新は一般には使えません。これは少なくとも \(f\) の勾配の Lipschitz 連続性、すなわち \(L\)-smoothness を仮定した近接勾配型の式であり、単に \(f\) が Lipschitz 連続であることからは導けません。

\(f,g\) が非平滑な Lipschitz 凸関数である場合、基本となる選択肢は

\[
\boxed{\text{厳密DCA} \quad+\quad \text{非平滑凸子問題ソルバ}}
\]

です。

## 1. 非平滑の場合の標準DCA

\(x_k\in X\) に対して
\[
y_k\in\partial g(x_k)
\]
を取得し、次の凸問題を解きます。

\[
x_{k+1}\in\arg\min_{x\in X}
\left\{
f(x)-\langle y_k,x\rangle
\right\}.
\tag{1}
\]

これは \(f\) の微分可能性を全く必要としません。\(f\) の値・劣勾配オラクルを用いて、凸非平滑最小化問題 (1) を解けばよいだけです。

凸性により
\[
g(x)\geq g(x_k)+\langle y_k,x-x_k\rangle
\]
なので、
\[
f(x)-g(x)
\leq
f(x)-g(x_k)-\langle y_k,x-x_k\rangle
\]
という意味で、\(g\) の線形化を用いたDCAが成立します。厳密に子問題を解けば、通常の条件の下で目的関数値は単調非増加です。DCAの標準形式は、まさに現在点で \(g\) の劣勾配を取り、得られた凸問題を解くものです。 [orsj](https://orsj.org/wp-content/corsj/or69-3/or69_3_136.pdf)

ただし、(1)を有限回の劣勾配評価だけで厳密に解けるとは限りません。したがって実際には、次の非厳密DCAが中心になります。

## 2. 非厳密DCA

子問題 (1) を反復法で近似的に解き、
\[
\widehat x_{k+1}\approx
\arg\min_{x\in X}
\{f(x)-\langle y_k,x\rangle\}
\]
とします。

自然な停止条件は、子問題
\[
\phi_k(x)=f(x)-\langle y_k,x\rangle+\iota_X(x)
\]
に対する最適性残差
\[
\operatorname{dist}\bigl(0,\partial\phi_k(\widehat x_{k+1})\bigr)
\leq \eta_k
\tag{2}
\]
です。

ここで
\[
\partial\phi_k(x)
=
\partial f(x)-y_k+N_X(x)
\]
なので、\(f\) の劣勾配 \(s_{k+1}\in\partial f(\widehat x_{k+1})\) と、制約集合の法錐 \(v_{k+1}\in N_X(\widehat x_{k+1})\) を使って
\[
\|s_{k+1}-y_k+v_{k+1}\|\leq\eta_k
\]
を確認します。

典型的な誤差条件は
\[
\eta_k\to0
\qquad\text{または}\qquad
\sum_{k=0}^{\infty}\eta_k<+\infty.
\]
この場合、適切な有界性や閉性の条件の下で、極限点について
\[
0\in \partial f(\bar x)-\partial g(\bar x)+N_X(\bar x)
\]
というDC臨界性が得られます。

## 3. 子問題の現実的な解法

### 3.1 劣勾配法

子問題
\[
\min_{x\in X}\phi_k(x)
=
\min_{x\in X}\{f(x)-\langle y_k,x\rangle\}
\]
に対して、
\[
s_t\in\partial f(z_t)
\]
を取り、
\[
z_{t+1}
=
\Pi_X\bigl(z_t-\alpha_t(s_t-y_k)\bigr)
\tag{3}
\]
とします。

ここで \(f\) が \(L_f\)-Lipschitz なら、通常
\[
\|s_t\|\leq L_f
\]
と評価できます。ただし、子問題の劣勾配は
\[
s_t-y_k
\]
なので、\(y_k\) も有界である必要があります。例えば \(g\) が \(L_g\)-Lipschitz なら
\[
\|y_k\|\leq L_g
\]
です。

したがって、
\[
\|s_t-y_k\|\leq L_f+L_g
\]
を使えます。標準的なステップ幅
\[
\alpha_t\asymp \frac{1}{\sqrt{t}}
\]
では、凸子問題の関数値誤差は一般に
\[
O(t^{-1/2})
\]
程度です。劣勾配法は汎用的ですが、子問題を高精度に解く用途では遅いので、通常は粗い初期段階の内側反復に使うのが現実的です。非平滑凸最小化に劣勾配法を使う基本形式は、例えば  に整理されています。 [msi.co](https://www.msi.co.jp/solution/nuopt/docs/glossary/articles/Subgradient.html)

### 3.2 Bundle法

\(f\) の値と劣勾配を利用できるなら、子問題にはbundle法が非常に自然です。

過去の点 \(z_i\) と劣勾配 \(s_i\in\partial f(z_i)\) から
\[
f(x)\geq f(z_i)+\langle s_i,x-z_i\rangle
\]
なので、下モデル
\[
\underline f_t(x)
=
\max_{i\in\mathcal B_t}
\left\{
f(z_i)+\langle s_i,x-z_i\rangle
\right\}
\]
を構成します。

子問題の近似として
\[
\min_{x\in X}
\left\{
\underline f_t(x)-\langle y_k,x\rangle
+\frac{\rho}{2}\|x-\bar z_t\|^2
\right\}
\tag{4}
\]
を解きます。

これは二次計画問題または線形計画問題に帰着しやすく、劣勾配法よりも実用的に高速なことがあります。特に、\(f\) が折れ線凸関数、最大アフィン関数、hinge loss、絶対値和などの場合に適しています。

DCAの外側反復とbundle法の内側反復を組み合わせるなら、

\[
\text{outer DCA}
\quad\longleftrightarrow\quad
\text{inner bundle solve}
\]

という二重ループになります。

### 3.3 ADMM・分解法

\(f\) が
\[
f(x)=f_1(Ax)+f_2(x)
\]
や
\[
f(x)=\sum_{j=1}^m f_j(x_j)
\]
のように分解できる場合、子問題 (1) を補助変数で分割してADMMを使えます。

例えば
\[
\min_{x,z}\;
f_1(z)+f_2(x)-\langle y_k,x\rangle
\quad
\text{s.t. } z=Ax,\ x\in X
\]
とし、各ブロックの非平滑凸問題を解きます。

この場合も \(C^2\) 性は不要です。必要なのは、各ブロックの近接問題または低次元凸問題が効率的に解けることです。

## 4. 近接正則化は使えるが、意味が違う

\(C^2\) 性がなくても、次のような近接正則化DCAは使えます。

\[
x_{k+1}\approx
\arg\min_{x\in X}
\left\{
f(x)-\langle y_k,x-x_k\rangle
+\frac{\rho_k}{2}\|x-x_k\|^2
\right\}.
\tag{5}
\]

ここで重要なのは、(5)を

\[
x_{k+1}
=
\Pi_X\left(
x_k-\frac{1}{\rho_k}
(\nabla f(x_k)-y_k)
\right)
\]

とは書けないことです。\(f\) は非平滑なので、(5)自体を**非平滑強凸凸最適化問題として解く**必要があります。

例えばその子問題に対し、劣勾配法を用いるなら
\[
s_t\in\partial f(z_t),
\]
\[
z_{t+1}
=
\Pi_X\left(
z_t-\alpha_t
\left[
s_t-y_k+\rho_k(z_t-x_k)
\right]
\right)
\tag{6}
\]
とします。

二次項のおかげで子問題は強凸になります。そのため、単なる子問題より安定しやすく、内側反復の停止基準も設定しやすいです。ただし、強凸でも \(f\) が非平滑なら、通常の勾配法のような \(O(1/t)\) の関数値収束を自動的に得られるわけではありません。bundle法、加速近接法、半滑らか性などの追加構造が必要です。

## 5. Bregman型について

前回挙げたBregman proximal DCAや外挿付きBregman法は、一般には

- 片方の凸関数が特定の滑らかさを持つ
- Bregman生成関数が適切に選べる
- 子問題が解ける

といった追加条件を仮定します。

したがって、今回のように \(f,g\) について「Lipschitz連続まで」しか仮定しないなら、標準的な \(L\)-smooth型のpDCAeを無条件に適用するのは適切ではありません。Bregman型も、単なるLipschitz性だけから直ちに利用できるわけではありません。 [link.springer](https://link.springer.com/article/10.1007/s10589-022-00411-w)

## 6. Lipschitz仮定から得られること

\(f,g\) が \(X\) 上でそれぞれ \(L_f,L_g\)-Lipschitzなら、劣勾配について概ね
\[
\|s\|\leq L_f,\qquad s\in\partial f(x),
\]
\[
\|y\|\leq L_g,\qquad y\in\partial g(x)
\]
が得られます。ただし、これは \(x\) が関数の定義域の相対内部などにあることを含めた標準的な条件の下で考えるべきです。

このため、DCAの線形化項は有界ですが、次の点には注意が必要です。

- \(X\) が非有界だと、子問題が非有界になる可能性がある。
- Lipschitz性だけでは、DCAの外側反復の一様な速度保証は一般には出ない。
- 子問題を厳密に解かない場合、近似誤差の管理が必要。
- DCAが得るのは一般には大域解ではなく、DC臨界点である。
- \(C^2\) 性がないため、Newton法や通常のNesterov型加速をそのまま使えない。

## 実際の推奨

今回の仮定なら、私は次の順序を推奨します。

### \(X\) への射影が安い場合

\[
\boxed{
\text{proximal DCA}
+
\text{内側のbundle法または劣勾配法}
}
\]

特に、子問題を完全に解かず、数回の内側反復で止める非厳密DCAが現実的です。

### \(f\) が最大アフィン・折れ線凸の場合

\[
\boxed{
\text{標準DCA}
+
\text{LP/QP型bundle子問題}
}
\]

この場合、劣勾配情報をbundleに蓄積できるため、単純な劣勾配法よりかなり有利です。

### 有限和構造がある場合

\[
\boxed{
\text{stochastic/incremental DCA}
}
\]

各 \(g_i\) または \(f_i\) の劣勾配をミニバッチで評価します。有限和構造を持つ非平滑DC問題に対する確率的DCAの研究もあります。 [jmlr](https://jmlr.org/papers/volume23/21-1146/21-1146.pdf)

### 構造が特にない場合

\[
\boxed{
\text{標準DCA}
+
\text{子問題に投影劣勾配法}
}
\]

ただし、速度は一般に遅く、理論上も主として
\[
0\in\partial f(\bar x)-\partial g(\bar x)+N_X(\bar x)
\]
という臨界性の保証になります。

したがって、今回の仮定での最も正確な表現は、

\[
\boxed{
y_k\in\partial g(x_k),\qquad
x_{k+1}\approx
\arg\min_{x\in X}
\{f(x)-\langle y_k,x\rangle\}
}
\]

であり、速度改善には \(C^2\) 性ではなく、

- \(X\) への射影可能性
- \(f\) の近接写像可能性
- 最大アフィン・分離可能・有限和などの構造
- bundleモデルの利用
- 子問題の強凸正則化

を利用する、という整理になります。
--->
