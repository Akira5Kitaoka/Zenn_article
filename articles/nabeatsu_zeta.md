---
title: "ナベアツゼータ関数入門"
emoji: "🔑"
type: "idea"  # tech or idea
topics: ["math", "数論", "複素解析"]
published: false
---

# これは何？

ナベアツ数から作られるゼータ関数の性質について調べた記事である．

## ノーテーション
- $\mathrm{NB}$をナベアツ数の集合，すなわち，3の倍数もしくは，3の付く数字（正整数）の集合とする．
- $\mathrm{NB}^c$を非ナベアツ数の集合，すなわち，3の倍数でも3がつくわけでもない数字（正整数）の集合とする．

# 定義とやりたいこと
ナベアツゼータ関数，非ナベアツゼータ関数を以下で定義する．

$$
    \zeta_{\mathrm{NB}} (s) := \sum_{n \in \mathrm{NB}} n^{-s},
    \quad 
    \zeta_{\mathrm{NB}^c} (s) := \sum_{n \in \mathrm{NB}^c} n^{-s}.
$$


お手本は母関数表示

$$
    \pi\cot\pi z=\frac1z+\sum_{n\ge1}\frac{2z}{z^2-n^2}=\frac1z-2\sum_{k\ge1}\zeta(2k)\,z^{2k-1}
$$

である．リーマンゼータ関数 $\zeta(s) := \sum_{n \geq 1} n^{-s}$ の偶数点の値 $\zeta(2k)$ を係数にもつ冪級数 $\sum_{k\ge1}\zeta(2k)\,z^{2k-1}$ は，閉じた関数 $\frac12\bigl(\frac1z-\pi\cot\pi z\bigr)$ で表される．同じことが**ナベアツゼータでもできるか？**——というのが本記事の問いである．

# 結果の要約

先に結論を述べる．

1. 3を桁に一つも含まない正整数の集合を $A$，$\omega=e^{2\pi i/3}$ とし，「3を桁に含まない数のゼータ」3つを

$$
F_j(s):=\sum_{n\in A}\omega^{jn}n^{-s}\qquad(j=0,1,2)
$$

とおくと，ナベアツゼータ・非ナベアツゼータは

$$
\zeta_{\mathrm{NB}^c}(s)=\frac23\bigl(F_0(s)-\operatorname{Re}F_1(s)\bigr),\qquad
\zeta_{\mathrm{NB}}(s)=\zeta(s)-\zeta_{\mathrm{NB}^c}(s)
$$

と，リーマンゼータ $\zeta$ と $F_j$ の明示的な結合で書ける（導出はステップ1）．

2. $\pi\cot\pi z$ の役割を果たす**非ナベアツ余接関数**が構成でき，非ナベアツゼータの偶数点の値はそのテイラー係数列として完全に特徴づけられる（構成はステップ2）．1. と合わせると，ナベアツゼータ関数も，リーマンゼータ関数と同様に母関数を持つ．


# 復習：リーマンゼータ関数の母関数表示

余接の部分分数展開

$$
\pi\cot\pi z=\frac1z+\sum_{n\ge1}\frac{2z}{z^2-n^2}=\frac1z-2\sum_{k\ge1}\zeta(2k)\,z^{2k-1}
$$

により，$\zeta(2k)$ たちは一つの関数 $\pi\cot\pi z$ のテイラー係数として束ねられる．(cf. [Jinbo2003, $\S$ 5.2])

本記事で行うのは，この**母関数による束ね上げのナベアツ世界への移植**である．以下，材料の分解から順に進める．

# ステップ1：ナベアツゼータを分解する

## 3を桁に含まない数のゼータ

3を桁に一つも含まない正整数の集合を $A$，使える桁の集合を $D=\{0,1,2,4,5,6,7,8,9\}$（9個）とし，

$$
F(s) := \sum_{n \in A} n^{-s}
$$

とおく．

## 包除原理

$M=\{n : 3\mid n\}$（3の倍数），$T=\{n : \text{3を桁に含む}\}$ とおくと $\mathrm{NB}=M\cup T$，$A=T^{c}$ である．$\zeta_3(s) := \sum_{3 \mid n} n^{-s} = 3^{-s}\zeta(s)$ とおくと，包除原理より

$$
\begin{align*}
\zeta_{\mathrm{NB}}(s)& =\zeta_3(s)+\bigl(\zeta(s)-F(s)\bigr)-\bigl(\zeta_3(s)-G(s)\bigr) \\
&=\zeta(s)-F(s)+G(s),
\\
\zeta_{\mathrm{NB}^c}(s)&=F(s)-G(s),
\end{align*}
$$

ここで

$$
G(s):=\sum_{\substack{n\in A\\ 3\mid n}} n^{-s}
$$

は「3を桁に含まない3の倍数」のゼータである．$\zeta_3$ が途中で相殺すること，両者の和が $\zeta$ に戻ることに注意する．

## $G$ を指標で分解する

残った $G$ は，$\bmod 3$ の加法指標で $F$ の仲間に分解できる．$\omega=e^{2\pi i/3}$ とし，

$$
F_j(s):=\sum_{n\in A}\omega^{jn}n^{-s}\qquad(j=0,1,2,\ F_0=F)
$$

とおくと，指標の直交性から

$$
G(s)=\frac13\bigl(F_0(s)+F_1(s)+F_2(s)\bigr),\qquad
\zeta_{\mathrm{NB}^c}(s)=\frac23\bigl(F_0(s)-\operatorname{Re}F_1(s)\bigr)
$$

（$F_2=\overline{F_1}$ を使った）．まとめると，目標のゼータたちは

$$
\zeta_{\mathrm{NB}^c}(s)=\frac23\bigl(F_0(s)-\operatorname{Re}F_1(s)\bigr),\qquad
\zeta_{\mathrm{NB}}(s)=\zeta(s)-\zeta_{\mathrm{NB}^c}(s)
$$

と，$\zeta$ と $F_j$ だけで書ける．あとは $F_j$ を計算できればよい．




# ステップ2：余接関数擬きの構成


## 練習
$\zeta(2k)$ たちを一つの関数のテイラー係数として束ね，その関数を関数等式で特徴づける——この仕組みは，$F$ に対しても再現できる．$A_0=\{0\}\cup A$ に対して

$$
C(z):=\sum_{n\in A_0}\frac{1}{z-n}
$$

と定義する．古典的な場合（$A_0=\mathbb{Z}_{\ge0}$）はこの和が発散するため $\pm n$ の対称化が必須だが，3を含まない数の集合では以下が成り立つ．

> **命題 1.** (1) $\sum_{n\in A}1/n  <\infty$．特に，(2) $C(z)$ を定める級数は $\mathbb{C} \setminus A_0$ で絶対収束し，そのコンパクト部分集合上で一様収束する．とくに $C$ は $\mathbb{C} \setminus A_0$ 上の正則関数を定める．

**証明**

(1) $L$ 桁の $A$ の元は $8\cdot9^{L-1}$ 個（先頭は $0,3$ 以外の8通り，残りは $3$ 以外の9通り）で，各元は $10^{L-1}$ 以上である．よって

$$
\sum_{n\in A}\frac1n\le\sum_{L\ge1}\frac{8\cdot9^{L-1}}{10^{L-1}}=\frac{8}{1-9/10}=80<\infty.
$$

(2) $z\in\mathbb{C}\setminus A_0$ を固定する．$n\ge2|z|$ なる $n\in A$ では $|z-n|\ge n-|z|\ge n/2$ ゆえ $\frac1{|z-n|}\le\frac2n$ であり，(1) より末尾の和は絶対収束する．残る有限個の項は有界なので，全体として絶対収束する．さらに $z$ が $\mathbb{C}\setminus A_0$ のコンパクト部分集合を動くとき，この評価は一様に取れるので級数は一様収束し，各項の正則性から $C$ の正則性が従う．Q.E.D.


以下の公式が成り立つ．
> **定理 1.**（乗法公式）　$C$ は
>
> $$
C(z)=\frac{1}{10}\sum_{d\in D}C\Bigl(\frac{z-d}{10}\Bigr)
$$
> 
> を満たす．

**証明**

各 $n\in A_0$ は末尾桁により $n=10m+d$（$m\in A_0,\ d\in D$）と一意に分解されるから，$A_0=\bigsqcup_{d\in D}(10A_0+d)$．よって

$$
\begin{align*}
C(z) & =\sum_{d\in D}\sum_{m\in A_0}\frac{1}{z-10m-d}
\\
& =\frac{1}{10}\sum_{d\in D}\sum_{m\in A_0}\frac{1}{\frac{z-d}{10}-m}\\
& =\frac{1}{10}\sum_{d\in D}C\Bigl(\frac{z-d}{10}\Bigr).
\end{align*}
$$

Q.E.D.

これは古典恒等式 $\pi\cot\pi z=\frac1{10}\sum_{d=0}^{9}\pi\cot\bigl(\pi\frac{z-d}{10}\bigr)$（余接の乗法公式）の，桁集合を $D$（3抜き）に取り替えたものに他ならない．しかも，この方程式が $C$ を特徴づける．

> **定理 2.**（一意性）　$\mathbb{C}$ 上の有理型関数で，極が $A_0$ の各点での留数 $1$ の単純極のみであり，かつ定理 1 の乗法公式を満たすものは，$C$ ただ一つである．

**証明**

2つの解の差を $h$ とすると，両者の極の位置と主要部（留数 $1$ の単純極）が一致しているので特異点はすべて除去可能であり，$h$ は整関数で，同じ方程式 $h(z)=\frac1{10}\sum_{d\in D}h\bigl(\frac{z-d}{10}\bigr)$ を満たす．

$M(R):=\max_{|z|\le R}|h(z)|$（**閉円板**上の最大値）とおく．$|z|\le R$ のとき，各 $d\in D$ に対して $\bigl|\frac{z-d}{10}\bigr|\le\frac{|z|+9}{10}\le\frac{R+9}{10}$ だから，方程式と三角不等式より

$$
|h(z)|\le\frac1{10}\sum_{d\in D}\Bigl|h\Bigl(\frac{z-d}{10}\Bigr)\Bigr|\le\frac{9}{10}\,M\Bigl(\frac{R+9}{10}\Bigr)
$$

（$D$ の元は $9$ 個しかないことに注意）．左辺の $|z|\le R$ にわたる最大値を取って

$$
M(R)\le\frac{9}{10}\,M\Bigl(\frac{R+9}{10}\Bigr).
$$

いま $R\ge1$ とすると $\frac{R+9}{10}\le R$ なので，$R_0=R,\ R_{j+1}=\frac{R_j+9}{10}$ で定めた列は $1\le R_j\le R$ を満たし続ける．上の不等式を $j$ 回反復し，$M$ の単調性（円板の包含 $R_j\le R$）を使うと

$$
M(R)\le\Bigl(\frac{9}{10}\Bigr)^{j}M(R_j)\le\Bigl(\frac{9}{10}\Bigr)^{j}M(R)\xrightarrow{\ j\to\infty\ }0.
$$

$M(R)$ は連続関数のコンパクト集合上の最大値ゆえ有限だから $M(R)=0$．$R\ge1$ は任意なので $h\equiv0$．Q.E.D.

なお，桁の個数が $9<10$ であることが縮小率 $\frac9{10}$ を与えている．全桁 $\{0,\dots,9\}$ の場合（$\pi\cot\pi z$ の世界）はこの因子が $1$ になってこの論法は使えず，実際に一意性も破れる：任意の定数 $c$ が斉次方程式 $c=\frac1{10}\sum_{d=0}^{9}c$ を満たすため，解に定数を足す自由度が残る．

そして，$C$ のテイラー係数がゼータの値そのものである．$|z|<1$ で $\frac{1}{z-n}=-\sum_{m\ge1}n^{-m}z^{m-1}$ と展開すれば（和の順序交換の正当化は後の命題 2 の証明と同様），

$$
C(z)=\frac1z-\sum_{m\ge1}F_0(m)\,z^{m-1},
$$

対称化すれば偶数点だけが残り，

$$
\frac{C(z)-C(-z)}{2}=\frac1z-\sum_{k\ge1}F_0(2k)\,z^{2k-1}
$$

となって，$\pi\cot\pi z=\frac1z-2\sum_{k\ge1}\zeta(2k)z^{2k-1}$ と完全に並行する．

## 非ナベアツ余接関数の定義

非ナベアツ版は $\bmod 3$ のクラス分けで得られる．関数を

$$
C_r(z):=\sum_{n\in A_0,\ n\equiv r\pmod3}(z-n)^{-1}
$$

とおく．$C_r(z)$ を定める級数は命題 1 と同様に $\mathbb{C} \setminus A_0$ 上で絶対収束し，

$$
    C (z) = C_0 (z) + C_1 (z) + C_2 (z)
$$

である．$10\equiv1\pmod3$ より連立乗法公式

$$
C_r(z)=\frac1{10}\sum_{d\in D}C_{\,r-d\ \mathrm{mod}\ 3}\Bigl(\frac{z-d}{10}\Bigr)\qquad(r=0,1,2)
$$

が成り立つ（証明は定理 1 と同じ）．注意として，右辺には $C_0$ も現れるため，スカラーの乗法公式（定理 1）を満たすのは $C$ 自身だけであり，個々の $C_r$ やその部分和は満たさない．一方，一意性は連立系のレベルで成り立つ：定理 2 の証明で $M(R):=\max_{r}\max_{|z|\le R}|h_r(z)|$ と取り替えれば同じ縮小論法が通り，極の条件＋連立乗法公式が三つ組 $(C_0,C_1,C_2)$ を一意に特徴づける．

この三つ組から $C^{\mathrm{NB}^c}:=C_1+C_2$ とおいたものが**非ナベアツ余接関数**である．実際，$C^{\mathrm{NB}^c}$ のテイラー係数が非ナベアツゼータの値を束ねている：

> **命題 2.**　$C^{\mathrm{NB}^c}$ は $|z|<1$ で正則であり，
>
> $$
\frac{C^{\mathrm{NB}^c}(z)-C^{\mathrm{NB}^c}(-z)}{2}=-\sum_{k\ge1}\zeta_{\mathrm{NB}^c}(2k)\,z^{2k-1}$$
> 
> が成り立つ．

**証明**

まず，$A_0$ の元のうち $\bmod 3$ でクラス $1,2$ に属するものは，ちょうど非ナベアツ数である．実際，$n\in A_0$ かつ $n\equiv1,2\pmod3$ なら，$n\neq0$，$n$ は3を桁に含まず，かつ $3\nmid n$，すなわち $n\in\mathrm{NB}^c$．逆に $n\in\mathrm{NB}^c$ なら $n\in A_0$ かつ $3\nmid n$．よって

$$
C^{\mathrm{NB}^c}(z)=C_1(z)+C_2(z)=\sum_{n\in\mathrm{NB}^c}\frac{1}{z-n}.
$$

$\min\mathrm{NB}^c=1$ なので（$0$ はクラス $0$ に属し，この和に現れない），右辺の各項は $|z|<1$ で極を持たず，$C^{\mathrm{NB}^c}$ はそこで正則である．$|z|<1$ で各項を等比級数に展開すると

$$
\frac{1}{z-n}=-\frac{1}{n}\cdot\frac{1}{1-z/n}=-\sum_{m\ge1}n^{-m}z^{m-1}
$$

であり，二重和は $n\ge1$ より

$$
\sum_{n\in\mathrm{NB}^c}\sum_{m\ge1}n^{-m}|z|^{m-1}
=\sum_{n\in\mathrm{NB}^c}\frac{1}{n}\cdot\frac{1}{1-|z|/n}
\le\frac{1}{1-|z|}\sum_{n\in\mathrm{NB}^c}\frac{1}{n}<\infty
$$

（最後は命題 1 (1) と $\mathrm{NB}^c\subset A$ による）と絶対収束するから，和の順序を交換してよく，

$$
C^{\mathrm{NB}^c}(z)=-\sum_{m\ge1}\Bigl(\sum_{n\in\mathrm{NB}^c}n^{-m}\Bigr)z^{m-1}
=-\sum_{m\ge1}\zeta_{\mathrm{NB}^c}(m)\,z^{m-1}.
$$

$z$ を $-z$ に置き換えると $C^{\mathrm{NB}^c}(-z)=-\sum_{m\ge1}\zeta_{\mathrm{NB}^c}(m)\,(-1)^{m-1}z^{m-1}$．差を取って $2$ で割ると，$m$ が奇数の項は打ち消し合い，$m=2k$（偶数）の項だけが残って

$$
\frac{C^{\mathrm{NB}^c}(z)-C^{\mathrm{NB}^c}(-z)}{2}=-\sum_{k\ge1}\zeta_{\mathrm{NB}^c}(2k)\,z^{2k-1}.
$$

Q.E.D.



# まとめ

オイラーの公式との対応を表にすると：

| | リーマンゼータ関数 | 非ナベアツゼータ関数 |
|---|---|---|
| 母関数 | $\dfrac12\Bigl(\dfrac1z-\pi\cot\pi z\Bigr)$ | $-\dfrac{C^{\mathrm{NB}^c}(z)-C^{\mathrm{NB}^c}(-z)}{2}$（$C^{\mathrm{NB}^c}=C_1+C_2$） |
| 特徴づけ | 周期性＋極 | 連立乗法公式＋極が $(C_0,C_1,C_2)$ を特徴づける |
| $z^{2k-1}$ の係数 | $\zeta(2k)=\frac{(-1)^{k+1}B_{2k}(2\pi)^{2k}}{2(2k)!}$ | $\zeta_{\mathrm{NB}^c}(2k)$ |

リーマンゼータを支配するのが周期性（$\bmod$ の世界）であるのに対し，ナベアツ数を支配するのは10進の桁の自己相似性である．アホになる数を数えていたら，自己相似性を持つ新しい特殊関数に辿り着いた，というのが本記事の結論である．

# 謝辞
Fable 5.

# 参考文献

[Jinbo2003] 神保道夫，複素解析入門，岩波書店，2003年．