anova.fr <- function(a, b, c, d){
  fac <- rep(rep(1:2), c(length(a), length(c)))
  ac <- c(a, c)
  bd <- c(b, d)
  data <- matrix(c(fac, ac, bd), length(ac), byrow=F)
  datafr <- data.frame(data)
  anovakun(datafr, "AsB", 2, 2)
}

#【ANOVA君：要因計画のタイプと水準数を入力することにより，分散分析を行う関数】
#1)フリーの統計ソフトウェア「R」で動作する関数
#2)被験者間要因（独立測度），被験者内要因（反復測度）のいずれか，または，両方を含む各タイプの分散分析を扱う
#3)引数としては，最初にデータフレーム名，次に計画のタイプ（""で囲むこと）を入力し，その後，各要因の水準数を順に入力する
#（作成：井関龍太）
#
#【使用法】
#anovakunに読み込むためのデータフレームは，以下のような形式で作っておく
#1)被験者間要因はタテ，被験者内要因はヨコに並べる
#2)被験者内要因を表す行はデータとして読み込まない（下の例では，点線より上の行は読み込まず，点線より下のデータのみ読み込む）
#3)被験者間要因を表す列はデータとして読み込む（下の例では，a1などの文字を含む列；被験者間要因数が増えるたびにラベル用の列を増やす）
#4)被験者間要因を表すラベルは，水準ごとに別の文字または数字を用いる（同じラベル＝同じ水準と見なされる）
#5)被験者内，被験者間要因ともに，前の方の要因から順に各水準のデータを入れ子状に整理して並べること（例を参照）
#6)下の被験者１～６の列は例の説明のためにつけたものなので，実際のデータフレームには必要ない
#
#［AsBC計画の例］（被験者間要因でデータ数が同じ）
#		b1	b1	b2	b2
#		c1	c2	c1	c2
#	------------------------------------
#	a1	12	9	14	13	---被験者１
#	a1	13	10	14	12	---被験者２
#	a1	11	10	13	15	---被験者３
#	a2	18	12	16	15	---被験者４
#	a2	17	14	15	14	---被験者５
#	a2	15	13	18	15	---被験者６
#
#データフレームを代入した変数名をxとすると，
#
#> anovakun(x, "AsBC", 2, 2, 2)
#
#のようにして関数を呼び出す
#
#［ABsC計画の例］（被験者間要因でデータ数がふぞろい）
#			c1	c2
#	-------------------------------
#	a1	b1	3.5	4.2	---被験者１
#	a1	b1	2.7	3.2	---被験者２
#	a1	b2	2.5	3.8	---被験者３
#	a1	b3	4.0	3.9	---被験者４
#	a2	b1	3.3	4.0	---被験者５
#	a2	b1	1.4	2.5	---被験者６
#	a2	b2	3.7	4.2	---被験者７
#	a2	b2	2.2	4.2	---被験者８
#	a2	b3	1.3	2.1	---被験者９
#	a2	b3	3.4	3.9	---被験者10
#
#データフレームを代入した変数名をxとすると，
#
#> anovakun(x, "ABsC", 2, 3, 2)
#
#のようにして関数を呼び出す
#
#【オプション】
#1)type2……type2 = Tとすると，平方和の計算方法をタイプⅡに切り替える
#2)tech……テクニカルアウトプット；tech = Tとすると，データフレームをリストでつないだ形式で結果を出力する
#3)data.frame……data.frame = Tとすると，計算に使用したデータフレームを出力する（関数中でdatと表現されているデータフレーム）
#4)holm……holm = Tとすると，多重比較の方法がHolmの方法になる
#5)hc……hc = Tとすると，多重比較の方法がHolland-Copenhaverの方法になる
#6)s2r……s2r = Tとすると，Shafferの方法のための仮説数の計算法を具体的な棄却のパターンに基づく方法に変更する（Rasmussenのアルゴリズムに基づく）
#7)s2d……s2d = Tとすると，Shafferの方法のための仮説数の計算法を具体的な棄却のパターンに基づく方法に変更する（Donoghueのアルゴリズムに基づく）
#8)fs1……fs1 = Tとすると，ステップ１の基準をステップ２の基準に置き換えた方法でShafferの方法を行う
#9)fs2r……fs2r = Tとすると，Shaffer2の方法とF-Shaffer1の方法を組み合わせた方法でShafferの方法を行う（Rasmussenのアルゴリズムに基づく）
#10)fs2d……fs2d = Tとすると，Shaffer2の方法とF-Shaffer1の方法を組み合わせた方法でShafferの方法を行う（Donoghueのアルゴリズムに基づく）
#11)hc, s2r……hc = Tかつs2r = T（または，fs2r = T，s2d = T，fs2d = T）とすると，Shaffer2の基準を用いてHolland-Copenhaverの方法を行う
#12)holm, hc……holm = Tかつhc = Tとすると，Holmの調整基準にSidakの不等式を用いた多重比較（Holm-Sidak法）を行う
#13)criteria……criteria = Tとすると，多重比較の出力において，調整済みｐ値の変わりに調整済みの有意水準を表示する
#14)lb……lb = Tとすると，すべての被験者内効果に対してイプシロンの下限値（Lower Bound）による自由度の調整を行う（Geisser-Greenhouseの保守的検定）
#15)gg……gg = Tとすると，すべての被験者内効果に対してGreehnouse-Geisserのイプシロンによる自由度の調整を行う
#16)hf……hf = Tとすると，すべての被験者内効果に対してHuynh-Feldtのイプシロンによる自由度の調整を行う
#17)auto……auto = Tとすると，球面性検定が有意であった被験者内効果に対してGreehnouse-Geisserのイプシロンによる自由度の調整を行う
#18)mau……mau = Tとすると，球面性検定の方法がMauchlyの球面性検定になる
#19)har……har = Tとすると，球面性検定の方法がHarrisの多標本球面性検定になる
#20)iga……iga = Tとすると，改良版一般近似検定を行う；イプシロンの代わりに各種の推定値を算出し，分散分析の際に適用する
#21)ciga……ciga = Tとすると，修正改良版一般近似検定を行う；イプシロンの代わりに各種の推定値を算出し，分散分析の際に適用する
#22)peta……peta = Tとすると，分散分析表に偏イータ二乗を追加する
#23)pomega……pomega = Tとすると，分散分析表に偏オメガ二乗を追加する
#24)prep……prep = Tとすると，分散分析表にp_repを追加する
#
#［オプション使用の例］（テクニカルアウトプットによる出力とHolmの方法による多重比較を指定）
#
#> anovakun(x, "AsB", 2, 2, tech = T, holm = T)
#
#【技術情報】
#1)ss.calc関数は，デフォルトでは，タイプⅢ平方和の計算法に基づいて分散分析を行う
#2)pro.fraction関数は，単純主効果の検定において誤差項をプールしない（水準別誤差項を使用；サブセットに分散分析を再適用するのと同じ）
#3)mod.Bon関数は，デフォルトでは，Shafferの方法による多重比較を行う（任意の棄却パターンにおける可能な真の帰無仮説の最大数に基づく方法）
#4)mod.Bon関数におけるShafferの方法，Holland-Copenhaverの方法の有意水準の計算は，Rasmussen（1993），Donoghue（2004）のアルゴリズムに基づく（オプションによる）
#5)mod.Bon関数による多重比較では，p値の低い対から順に基準値よりも値が低いか否かを判定し，いったん有意でない対が見つかったら，
#以降の対についてはすべて判断を保留する（p値が基準値を下回っても*マークを表示しない）；調整済みｐ値の表示の際には，調整済みｐ値が
#既定の有意水準（５％）を下回った対にのみ*マークを表示する
#6)epsilon.calc関数は，デフォルトでは，被験者内要因を含むデータに対してMendozaの多標本球面性検定を行う（近似カイ二乗による）
#7)epsilon.calc関数は，オプション指定により，被験者内要因を含むデータに対してMauchlyの球面性検定を行う（近似カイ二乗による）
#8)epsilon.calc関数は，オプション指定により，被験者内要因を含むデータに対してHarrisの多標本球面性検定を行う（近似カイ二乗による）
#9)anovakunを構成する関数は，仕様上は，最大で26要因までの計画に対応できる；この上限は，要因を表すラベルとしてアルファベット26文字（LETTERSとletters）を使用していることによる
#10)epsilon.calc関数によるIGAとCIGAは，非加重平均を想定した計算法に基づく（Algina, 1997; Algina & Oshima, 1995）
#11)anova.modeler関数は，IGAとCIGAを用いた際には，b_hat，c_hatによって調整した後のF値をF値の列に表示する
#
#【このファイルの含む関数】
#hmean……調和平均を計算する
#anovakun……データフレームの作成とプロセス全体の制御を行う
#elematch……文字列中のすべての要素を含む文字列をマッチングする
#ss.calc……平方和を計算する
#sig.sign……有意水準に合わせて記号を表示する
#epsilon.calc……Greenhouse-GeisserとHuynh-Feldtのイプシロンを計算する
#anova.modeler……分散分析表を作る
#mod.Bon……修正Bonferroniの方法による多重比較を行う
#post.analyses……下位検定を行う
#pro.fraction……効果のタイプに適した下位検定を割り当てる
#info.out……基本情報を出力する
#bstat.out……記述統計量を出力する
#anovatab.out……分散分析表を出力する
#mod.Bon.out……修正Bonferroniの方法による多重比較の結果を出力する
#simeffects.out……単純主効果の検定の結果を出力する
#post.out……下位検定の結果を出力する
#each.out……効果のタイプによって異なる出力方式を割り当てる
#
#【バージョン情報】
#1)anovakun version 1.0.0（R 2.5.1作成；2007/09/03公開）
#・分散分析，単純主効果の検定，多重比較
#2)anovakun version 2.0.0（R 2.5.1作成；2007/10/01公開）
#・球面性検定とイプシロン調整；epsilon.calc関数の追加とそれに伴う改修
#・平方和の計算方法を修正；ss.calc関数の修正とelematch関数の追加，それに伴う改修
#・多重比較を行う際に，他の要因の効果を考慮した上でのMSeを用いてｔ統計量を計算するように修正；mod.Bon関数と関連する部分の改修
#3)anovakun version 2.1.0（R 2.5.1作成；2007/11/01公開）
#・平方和の計算方法を変更し，高速化を試みる；ss.calc関数の変更；その他，最適化
#4)anovakun version 2.2.0（R 2.5.1，R 2.6.0作成；2007/12/03公開）
#・QR分解の方法をLAPACKに変更し，高速化を試みる；その他の修正
#5)anovakun version 3.0.0（R 2.5.1，R 2.6.0作成；2008/01/04公開）
#・タイプⅢ平方和を計算する機能を追加（タイプⅢをデフォルトに設定）；ss.calc関数の変更とそれに伴う改修
#・Shafferの方法のための可能な真の帰無仮説の数の計算アルゴリズムを追加；mod.Bon関数の変更とshaffer2，fshaffer1，fshaffer2オプションの追加
#・多重比較の出力において調整済みｐ値を表示する機能を追加（デフォルトに設定）；mod.Bon関数の変更とそれに伴う改修
#6)anovakun version 3.1.0（R 2.5.1，R 2.6.0作成；2008/02/01公開）
#・aov関数のアルゴリズムを利用して誤差平方和の計算の高速化を試みる；ss.calc関数の変更；その他の修正
#7)anovakun version 3.2.0（R 2.5.1，R 2.6.0作成；2008/04/01公開）
#・Shafferの方法のための別バージョンのアルゴリズムを追加；mod.Bon関数の変更とs2d，fs2dオプションの追加
#・オプション名の変更；“shaffer2，fshaffer1，fshaffer2”を“s2r，fs1，fs2r”に変更
#・効果量とp_repを計算する機能を追加；anova.modeler関数，anovatab.out関数ほかの変更
#8)anovakun version 4.0.0（R 2.5.1，R 2.6.0，R2.7.0作成；2008/06/02公開）
#・Mendozaの多標本球面性検定とHarrisの多標本球面性検定を行う機能を追加（Mendozaをデフォルトに設定）；epsilon.calc関数と関連する部分の改修
#・Huynhの改良版近似検定とAlgina-Lecoutreの修正改良版近似検定を行う機能を追加；epsilon.calc関数及び関連する箇所の改修
#・取り扱い可能な要因の数を拡張；anovakun関数，anova.modeler関数，epsilon.calc関数，pro.fraction関数ほかの改修
#・複数の効果量をオプション指定したときに，すべての指標が同時に出力されるように変更；anova.modeler関数，pro.fraction関数ほかの改修
#
#【この関数の使用に関して】
#1)anovakunとこれを構成する関数（コード群）は，自由に使用，改変，再配布していただいて結構です。
#2)ただし，改変を加えたものを公開する際には，改変を加えたことを明記し，メインの関数の名前をanovakun以外の名前に変更してください。
#3)anovakunとこれを構成する関数（コード群）の使用によって生じるいかなる結果に関しても作成者は責任を負いかねますのでご了承ください。


#調和平均を計算する関数
hmean <- function(datvector){
	length(datvector)/sum(1/datvector)
}


#ANOVA君本体：データフレームの作成，プロセス全体の制御を行う関数
anovakun <- function(dataset, design, ..., type2 = FALSE, tech = FALSE, data.frame = FALSE, holm = FALSE, hc = FALSE, 
	s2r = FALSE, s2d = FALSE, fs1 = FALSE, fs2r = FALSE, fs2d = FALSE, criteria = FALSE, lb = FALSE, gg = FALSE, hf = FALSE, 
	auto = FALSE, mau = FALSE, har = FALSE, iga = FALSE, ciga = FALSE, peta = FALSE, pomega = FALSE, prep = FALSE) {
	bet.with <- strsplit(design, "s")
	betN <- nchar(bet.with[[1]])[1]#被験者間要因がないときは“０”
	withN <- nchar(bet.with[[1]])[2]#被験者内要因がないときは“NA”
	maxfact <- nchar(design) - 1#実験計画全体における要因数
	eachlev <- unlist(list(...))#各要因の水準数

	if(betN == 0){#被験者内計画の場合
		if(withN == 1){#１要因の場合
			eachcue <- 1
			timescue <- 1
		}else{#その他の場合
			eachcue <- c(sapply(2:withN, function(x) prod(eachlev[x:withN])), 1)
			timescue <- prod(eachlev) / (eachcue * eachlev)
		}

		#各被験者内要因の各水準を表すベクトル
		mainlab <- as.data.frame(sapply(1:withN, function(x) rep(paste(letters[x], 1:eachlev[x], sep = ""), 
			each = nrow(dataset) * eachcue[x], times = timescue[x])))

	}else if(is.na(withN)){#被験者間計画の場合
		eachcue <- lapply(1:betN, function(x) as.vector(table(dataset[,x:1])))
		timescue <- c(1, sapply(1:(betN-1), function(x) prod(eachlev[1:x])))

		#各被験者間要因の各水準を表すベクトル
		mainlab <- as.data.frame(sapply(1:betN, function(x) rep(rep(paste(letters[x], 1:eachlev[x], sep = ""), 
			times = timescue[x]), times = eachcue[[x]])))

	}else{#混合要因計画の場合
		eachcue1 <- lapply(1:betN, function(x) as.vector(table(dataset[,x:1])))
		timescue1 <- c(1, sapply(1:(betN-1), function(x) prod(eachlev[1:x])))
		withlev <- prod(eachlev[-(1:betN)])#被験者内要因のすべての水準数をかけたもの

		#各被験者間要因の各水準を表すベクトル
		mainlab1 <- as.data.frame(sapply(1:betN, function(x) rep(rep(rep(paste(letters[x], 1:eachlev[x], sep = ""), 
			times = timescue1[x]), times = eachcue1[[x]]), withlev)))

		if(withN == 1){#１要因の場合
			eachcue2 <- 1
			timescue2 <- 1
		}else{#その他の場合
			eachcue2 <- c(sapply(2:withN, function(x) prod(eachlev[(x+betN):(withN+betN)])), 1)
			timescue2 <- withlev / (eachcue2 * eachlev[-(1:betN)])
		}

		#各被験者内要因の各水準を表すベクトル
		mainlab2 <- as.data.frame(sapply(1:withN, function(x) rep(paste(letters[x+betN], 1:eachlev[x+betN], sep = ""), 
			each = nrow(dataset) * eachcue2[x], times = timescue2[x])))

		mainlab <- cbind(mainlab1, mainlab2)
	}

	names(mainlab) <- LETTERS[1:maxfact]#要因のラベルを列名に設定

	depv <- unlist(dataset[,(betN+1):ncol(dataset)])#従属変数のベクトル
	names(depv) <- NULL#もとの変数名を消去

	#データフレームにまとめる
	dat <- data.frame("s" = rep(paste("s", 1:nrow(dataset), sep = ""), ncol(dataset) - betN), mainlab, "y" = depv)

	#作成したデータフレームをもとに，各条件ごとの平均と標準偏差を計算する
	sncol <- as.vector(tapply(dat$y, dat[,(maxfact+1):2], length))#セルごとのデータ数を計算
	mncol <- as.vector(tapply(dat$y, dat[,(maxfact+1):2], mean))#セルごとの平均を計算
	sdcol <- as.vector(tapply(dat$y, dat[,(maxfact+1):2], sd))#セルごとの標準偏差を計算

	#記述統計量の表において各要因の各水準を表すためのラベル（数字列）を作成
	maincols <- expand.grid(lapply((maxfact+1):2, function(x) levels(dat[,x])))
	maincols <- maincols[,order(maxfact:1)]#アルファベット順に並べ替え

	#記述統計量をデータフレームにまとめる
	bstatist <- data.frame(maincols, sncol, mncol, sdcol)
	names(bstatist) <- c(LETTERS[1:maxfact], "N", "Mean", "S.D.")#要因のラベルほかを列名に設定

	#記述統計量の表を出力する際の改行位置を指定する
	if(maxfact < 3) margin <- prod(eachlev)
	else margin <- prod(eachlev[-(1:(maxfact-2))])

	#anova.modelerにデータフレームを送り，分散分析の結果を得る
	mainresults <- anova.modeler(dat = dat, design = design, type2 = type2, 
		lb = lb, gg = gg, hf = hf, auto = auto, mau = mau, har = har, iga = iga, ciga = ciga, 
		peta = peta, pomega = pomega, prep = prep)

	#post.analysesにデータフレームと分散分析の結果を送り，下位検定の結果を得る
	postresults <- post.analyses(dat = dat, design = design, mainresults = mainresults, type2 = type2, 
		holm = holm, hc = hc, s2r = s2r, s2d = s2d, fs1 = fs1, fs2r = fs2r, fs2d = fs2d, criteria = criteria, 
		lb = lb, gg = gg, hf = hf, auto = auto, mau = mau, har = har, iga = iga, ciga = ciga, 
		peta = peta, pomega = pomega, prep = prep)

	#基本情報の取得
	info1 <- paste("[ ", design, "-Type Design ]", sep = "")#要因計画の型
	info2 <- paste("This output was generated via anovakun 4.0.0 at ", substr(R.version$version.string, 1, 15), ".", sep = "")#バージョン情報など
	info3 <- paste("It was executed on ", date(), ".", sep = "")#実行日時

	#Unbalancedデザイン（データ数ふぞろい）の場合，プロンプトを追加
	if(length(unique(sncol)) != 1){
		if(sum(is.na(mainresults$ano.info1)) == 1){
			if(type2 == TRUE) mainresults$ano.info1 <- c("== This data is UNBALANCED!! ==", "== Type II SS is applied. ==")
			else mainresults$ano.info1 <- c("== This data is UNBALANCED!! ==", "== Type III SS is applied. ==")
		}else{
			if(type2 == TRUE) mainresults$ano.info1 <- append(mainresults$ano.info1, c("== This data is UNBALANCED!! ==", "== Type II SS is applied. =="))
			else mainresults$ano.info1 <- append(mainresults$ano.info1, c("== This data is UNBALANCED!! ==", "== Type III SS is applied. =="))
		}
	}

	#結果を表示する
	if(tech == TRUE){#データフレーム形式での出力の場合
		if(data.frame == TRUE){return(list("INFORMATION" = rbind(info1, info2, info3), 
			"DESCRIPTIVE STATISTICS" = bstatist, 
			"SPHERICITY INDICES" = list(mainresults$epsi.info1, mainresults$epsitab), 
			"ANOVA TABLE" = list(mainresults$ano.info1, mainresults$anovatab), 
			"POST ANALYSES" = postresults,
			"DATA.FRAME" = dat))#計算に使用したデータフレームを付加
		}else{return(list("INFORMATION" = rbind(info1, info2, info3), 
			"DESCRIPTIVE STATISTICS" = bstatist, 
			"SPHERICITY INDICES" = list(mainresults$epsi.info1, mainresults$epsitab), 
			"ANOVA TABLE" = list(mainresults$ano.info1, mainresults$anovatab), 
			"POST ANALYSES" = postresults))
		}
	}else{#表形式での出力の場合
		if(data.frame == TRUE){info.out(info1, info2, info3)
			bstat.out(bstatist, maxfact, margin)
			anovatab.out(mainresults)
			post.out(postresults, design)
			list("DATA.FRAME" = dat)#計算に使用したデータフレームを付加
		}else{info.out(info1, info2, info3)
			bstat.out(bstatist, maxfact, margin)
			anovatab.out(mainresults)
			post.out(postresults, design)
		}
	}
}


#文字列中のすべての要素を含む文字列をマッチングする関数
#grepとの違いは，“A:C”などの文字列を照合パターンとした場合に，“A:B:C”のように間に別の文字を挟んだ文字列もマッチと判定する点
#照合パターンが１文字の場合は，grepと同じ結果を返す
elematch <- function(Mstrings, stex){
	#マッチングする文字列を分解して，それぞれgrep関数を適用
	matchlist <- lapply(strsplit(Mstrings, "")[[1]], function(x) grep(x, stex))

	#文字列の各要素とマッチした値の共通部分のみ取り出す
	buffer <- matchlist[[1]]
	if(length(matchlist) != 1){
		for(i in 2:length(matchlist)){
			buffer <- buffer[is.element(buffer, matchlist[[i]])]
		}
	}#文字列が１文字のときはgrepの結果をそのまま返す

	return(buffer)
}


#平方和を計算する関数
ss.calc <- function(full.elem, dat, type2 = FALSE){
	#full.elemの並べ替え
	elem.num <- seq(1, max(nchar(full.elem)), by = 2)
	full.elem <- unlist(sapply(elem.num, function(x) full.elem[nchar(full.elem) == x]))

	#計画行列を作る
	if(sum(grep("s", full.elem)) == 0) eff.elem <- full.elem
	else eff.elem <- full.elem[-grep("s", full.elem)]#誤差項との交互作用効果を除いたもの

	er.elem <- grep("s", full.elem, value = TRUE)#誤差項のみ
	eff.modeleq <- paste("~ ", gsub(",", " +", toString(eff.elem)), sep = "")
	er.modeleq <- paste("~ ", gsub(",", " +", toString(er.elem)), sep = "")

	def.contr <- options()$contrasts#contrastsのデフォルト設定を保存
	options(contrasts = c("contr.sum", "contr.sum"))#設定を変更
	dmat <- model.matrix(as.formula(eff.modeleq), dat)

	#計画行列とデータを統合する
	exmat <- as.matrix(cbind(dmat, dat$y))#拡大行列を作る
	promat <- crossprod(exmat)#拡大行列の積和行列
	endline <- nrow(promat)#積和行列の列数

	#各効果に対応する計画行列の列（行）番号を得る
	sepcol <- attr(dmat, "assign")#計画行列からピボットを表すベクトルを取り出す
	pivot.col <- lapply(1:max(sepcol), function(x) (1:length(sepcol))[sepcol == x])#各効果を表現する列の番号
	names(pivot.col) <- eff.elem

	if(type2 == TRUE){#線形モデルを用いてタイプⅡ平方和を計算する
		#各効果の平方和の計算
		#各モデルのための部分行列を選択
		ss.line1 <- lapply(eff.elem, function(x) c(1, unlist(pivot.col[match(x, names(pivot.col))]), unlist(pivot.col[-elematch(x, names(pivot.col))])))
		ss.line2 <- lapply(eff.elem, function(x) c(1, unlist(pivot.col[-elematch(x, names(pivot.col))])))

		#各モデルの最小二乗解を得てもとのベクトルにかけたものの和を取る
		ss.base1 <- sapply(ss.line1, function(x) sum(qr.coef(qr(promat[x, x], LAPACK = TRUE), promat[x, endline]) * promat[x, endline]))
		ss.base2 <- sapply(ss.line2, function(x) sum(qr.coef(qr(promat[x, x], LAPACK = TRUE), promat[x, endline]) * promat[x, endline]))

		#各効果を含むモデルと含まないモデルの差分を取る
		ss.all <- ss.base1 - ss.base2
		ss.all <- lapply(ss.all, function(x) x)#リスト化
		names(ss.all) <- eff.elem

	}else{#線形モデルを用いてタイプⅢ平方和を計算する
		#各効果の平方和の計算
		eff.line <- c(1, unlist(pivot.col))

		#各モデルのための部分行列を選択
		ss.line <- lapply(eff.elem, function(x) eff.line[-match(pivot.col[[match(x, names(pivot.col))]], eff.line)])

		#各モデルの最小二乗解を得てもとのベクトルにかけたものを合計
		ss.eff <- sum(qr.coef(qr(promat[eff.line, eff.line], LAPACK = TRUE), promat[eff.line, endline]) * promat[eff.line, endline])
		ss.base <- sapply(ss.line, function(x) sum(qr.coef(qr(promat[x, x], LAPACK = TRUE), promat[x, endline]) * promat[x, endline]))

		#各効果を含むモデルと含まないモデルの差分を取る
		ss.all <- ss.eff - ss.base
		ss.all <- lapply(ss.all, function(x) x)#リスト化
		names(ss.all) <- eff.elem
	}

	#全体平方和の計算
	ss.T <- promat[endline, endline] - sum(qr.coef(qr(promat[1, 1], LAPACK = TRUE), promat[1, endline]) * promat[1, endline])

	#誤差平方和の計算
	if(sum(grep("s", full.elem)) == 0){#誤差項を分解する必要がない場合
		ss.Er <- promat[endline, endline] - sum(qr.coef(qr(promat[1:(endline-1), 1:(endline-1)], LAPACK = TRUE), promat[1:(endline-1), endline]) * promat[1:(endline-1), endline])
	}else{#誤差項を分解する必要がある場合
		ss.Er <- NA
		emat <- model.matrix(as.formula(er.modeleq), dat)#切片と誤差項のみの計画行列
		er.num <- length(er.elem)#誤差項モデルの項の数

		qr.er <- qr(emat, LAPACK = TRUE)#誤差項モデルのQR分解
		er.col <- attr(qr.er$qr, "assign")[qr.er$pivot[1:qr.er$rank]]#有効ピボット

		tq.er <- t(qr.Q(qr.er))
		qty <- as.matrix(tq.er %*% dat$y)
		qtx <- tq.er %*% dmat

		ss.er <- rep(list(NA), er.num)#誤差平方和格納用のリストを宣言
		names(ss.er) <- er.elem

		for(i in 1:er.num){
			select <- er.col == i
			xi <- qtx[select, , drop = FALSE]
			cols <- colSums(xi^2) > 1e-05
			ss.er[[i]] <- sum(qr.resid(qr(xi[, cols, drop = FALSE]), qty[select, , drop = FALSE]) * qty[select, , drop = FALSE])
		}
		ss.all <- c(ss.all, ss.er)
	}

	ss.results <- c(ss.all, "ss.T" = ss.T, "ss.Er" = ss.Er)
	options(contrasts = def.contr)#contrastsの設定をもどす

	return(ss.results)

}


#有意水準に合わせて記号を表示する関数
sig.sign <- function(pvalue){
	ifelse(is.na(pvalue), "", 
	ifelse(pvalue < 0.001, "***", 
	ifelse(pvalue < 0.01, "**", 
	ifelse(pvalue < 0.05, "*", 
	ifelse(pvalue < 0.10, "+", "ns")))))
}


#Greenhouse-GeisserとHuynh-Feldtのイプシロンを計算する関数
#被験者内要因を含まない計画を投入すると適切に動作しないので注意
epsilon.calc <- function(dat, design, mau = FALSE, har = FALSE, iga = FALSE, ciga = FALSE){
	#要因計画の型から被験者内要因を特定
	bet.with <- strsplit(design, "s")

	#被験者間要因の特定
	othlabel <- match(strsplit(bet.with[[1]][1], "")[[1]], names(dat))
	othnum <- othlabel[order(othlabel, decreasing = TRUE)]#後の方の要因のラベルを先に並べる
	othmat <- sapply(othnum, function(x) nlevels(dat[,x]))#各要因の水準数をベクトル化
	ol <- ifelse(length(othmat) == 0, 1, prod(othmat))#全被験者間要因の組み合わせ水準数を取得；被験者間要因がなければ，１を代入

	#被験者内要因の特定
	replabel <- match(strsplit(bet.with[[1]][2], "")[[1]], names(dat))
	repnum <- replabel[order(replabel, decreasing = TRUE)]#後の方の要因のラベルを先に並べる
	repmat <- sapply(repnum, function(x) nlevels(dat[,x]))#各要因の水準数をベクトル化
	rl <- prod(repmat)#全被験者内要因の組み合わせ水準数を取得

	cellN <- length(unique(dat$s))#被験者間要因をつぶしての全被験者の数を取得

	if(length(othmat) == 0) othN <- cellN#被験者間要因がないときはcellNと同じ値を代入
	else othN <- as.vector(table(dat[names(dat)[othnum]]) / rl)#被験者間要因の各組み合わせにおける被験者数をベクトル化

	#データフレームを分割し，共分散行列を作る
	if(length(othnum) == 0){#被験者間要因がないときはデータフレームを分割しない
		covmatrices <- cov(as.data.frame(split(dat$y, dat[names(dat)[repnum]])))
	}else{#被験者間要因の組み合わせ水準ごとにデータフレームを分割
		covmatrices <- lapply(split(dat, dat[names(dat)[othnum]]), function(x) cov(as.data.frame(split(x$y, x[names(x)[repnum]]))))
	}

	#データフレームのリストを三次元配列に変換
	covbuffer <- array(unlist(covmatrices), dim = c(rl, rl, ol))

	#複数の共分散行列をプール
	tm <- apply(covbuffer, c(1, 2), function(x) sum((othN-1) * x)) / (cellN - ol)

	#正規直交対比行列を作る；被験者内要因の数に合わせて異なるパターンを得る
	replev <- sapply(replabel, function(x) nlevels(dat[,x]))#計画タイプの順に各被験者内要因の水準数を得る
	repleng <- length(replabel)#反復測定要因の数

	def.contr <- options()$contrasts#contrastsのデフォルト設定を保存
	options(contrasts = c("contr.helmert", "contr.helmert"))#設定を変更

	ortho.model <- expand.grid(lapply(repnum, function(x) levels(dat[,x])))
	names(ortho.model) <- paste("R", repleng:1, sep = "")
	ortho.helm <- model.matrix(as.formula(paste("~ ", gsub(",", " *", toString(paste("R", 1:repleng, sep = ""))), sep = "")), ortho.model)
	matcue <- attr(ortho.helm, "assign")
	matdivider <- sapply(1:max(matcue), function(x) length(matcue[matcue == x]))#正規直交行列を分割する際の行数

	options(contrasts = def.contr)#設定をもどす
	effect.name <- unlist(sapply(1:repleng, function(y) combn(names(dat)[replabel], y, function(x) gsub(", ", "x", toString(x)))))#効果のラベルを作る

	ortho.coef <- t(ortho.helm[, 2:rl, drop = FALSE])#直交対比のパターンのみを取り出す；行が１のときにベクトルに変換されないようにdrop = FALSEを使う
	ortho.denomi <- rowSums(ortho.coef^2)^(1/2)

	#パターンを直交対比行列にする
	orthoM <- ortho.coef / ortho.denomi

	#共分散行列と正規直交対比行列をかける
	otoM <- orthoM %*% tm %*% t(orthoM)

	#行列全体を使っての球面性検定
	if(mau == TRUE){
		#プロンプトの準備
		epsi.info1 <- paste("== Mauchly's Sphericity Test and Epsilons ==", sep = "")
		lamlab <- "W"

		#Mauchlyの球面性検定
		eps.Lambda <- det(otoM) / (sum(diag(otoM)) / (rl - 1))^(rl - 1)
		eps.m <- 1 - (2 * rl^2 - 3 * rl + 3) / (6 * (cellN - ol) * (rl - 1))

		if(cellN <= rl){#被験者数が被験者内要因の組み合わせ水準数を下回るときは妥当なカイ二乗値を計算できない
			epsChi <- NA
			epsi.info1 <- paste(epsi.info1, "\n", 
				"*** CAUTION! The test of GLOBAL SPHERICITY is INVALID because of small sample size. ***", "\n", 
				"*** The minimum sample size for valid computation is N = ", rl + 1, " at each group. ***", sep = "")
		}else{
			epsChi <- -(cellN - ol) * eps.m * log(eps.Lambda)
		}

		eps.df <- (((rl - 1) * rl) / 2) - 1
		eps.p <- pchisq(epsChi, ifelse(eps.df == 0, NA, eps.df), lower.tail = FALSE)
	}else if(har == TRUE){
		#プロンプトの準備
		epsi.info1 <- paste("== Harris's Multisample Sphericity Test and Epsilons ==", sep = "")
		lamlab <- "h_hat"

		#Harrisの多標本球面性検定
		proA <- array(apply(covbuffer, 3, function(x) orthoM %*% x %*% t(orthoM)), dim = c(rl-1, rl-1, ol))

		if(cellN <= rl){#被験者数が被験者内要因の組み合わせ水準数を下回るときは妥当なカイ二乗値を計算できない
			epsChi <- NA
			eps.Lambda <- NA
			epsi.info1 <- paste(epsi.info1, "\n", 
				"*** CAUTION! The test of GLOBAL SPHERICITY is INVALID because of small sample size. ***", "\n", 
				"*** The minimum sample size for valid computation is N = ", rl + 1, " at each group. ***", sep = "")
		}else{
			harTr <- apply(proA, 3, function(x) sum(diag(x)))
			harSq <- apply(proA, 3, function(x) sum(diag(x %*% x)))
			eps.Lambda <- sum((othN - 1) * harTr)^2 / sum((othN - 1) * harSq)
			epsChi <- pmax(0, ((cellN - ol) * (rl - 1) / 2) * ((cellN - ol) * (rl - 1) / eps.Lambda - 1))#負の値は０にそろえる
			eps.df <- (ol * (rl - 1) * rl) / 2 - 1
			eps.p <- pchisq(epsChi, ifelse(eps.df == 0, NA, eps.df), lower.tail = FALSE)
		}
	}else{
		#プロンプトの準備
		epsi.info1 <- paste("== Mendoza's Multisample Sphericity Test and Epsilons ==", sep = "")
		lamlab <- "Lambda"

		#Mendozaの多標本球面性検定
		ptm <- covbuffer * rep(othN, each = rl^2)
		proA <- array(apply(ptm, 3, function(x) orthoM %*% x %*% t(orthoM)), dim = c(rl-1, rl-1, ol))

		if(cellN <= rl){#被験者数が被験者内要因の組み合わせ水準数を下回るときは妥当なカイ二乗値を計算できない
			epsChi <- NA
			eps.Lambda <- NA
			epsi.info1 <- paste(epsi.info1, "\n", 
				"*** CAUTION! The test of GLOBAL SPHERICITY is INVALID because of small sample size. ***", "\n", 
				"*** The minimum sample size for valid computation is N = ", rl + 1, " at each group. ***", sep = "")
		}else{
			menL1 <- log((cellN-ol)^((cellN-ol)*(rl-1)/2) / prod((othN-1)^((othN-1)*(rl-1)/2)))
			menL2 <- apply(proA, 3, function(x) det(x))
			menL2 <- sum(log(ifelse(menL2 < 0, NA, menL2)) * (othN-1)/2)#行列式が負になったときは対数にできないのでNAを代入
			menL3 <- log(sum(diag(apply(proA, c(1, 2), function(x) sum(x))/(rl-1)))) * ((cellN-ol)*(rl-1)/2)
			menL <- menL1 + menL2 - menL3
			eps.m <- 1 - ((((cellN-ol) * (rl-1)^2 * rl * (2*(rl-1)+1) - (2*(cellN-ol)*(rl-1)^2)) * sum(1/(othN-1)) - 4) / (6*(cellN-ol)*(rl-1) * (ol*(rl-1)*rl-2)))
			eps.m[is.nan(eps.m)] <- 0#NaNが出たところには０を代入
			epsChi <- - 2 * eps.m * menL
			eps.Lambda <- exp(menL)
		}

		eps.df <- (ol * (rl - 1) * rl) / 2 - 1
		eps.p <- pchisq(epsChi, ifelse(eps.df == 0, NA, eps.df), lower.tail = FALSE)
	}

	#イプシロンを計算する
	LB.ep <- 1 / (rl - 1)
	GG.ep <- sum(diag(otoM))^2 / (nrow(otoM) * sum(otoM^2))
	HF.ep <- (cellN * (rl - 1) * GG.ep - 2) / ((rl - 1) * (cellN - ol - (rl - 1) * GG.ep))

	#被験者内要因の数によって処理を変更する
	if(repleng == 1){#被験者内要因が１つのときは有意性判定のマークを用意するのみ
		sig.mark <- sig.sign(eps.p)
		seportM <- list(orthoM)
	}else{#被験者内要因が複数あるときは各要因ごとに，検定統計量ととイプシロンを計算
		#ラベルの追加
		effect.name <- c("Global", effect.name)

		#正規直交対比行列を被験者内要因の水準数によって分割
		divpoint <- rbind(cumsum(matdivider) - (matdivider - 1), cumsum(matdivider))
		seportM <- unlist(apply(divpoint, 2, function(x) list(orthoM[x[1]:x[2], , drop = FALSE])), recursive = FALSE)
		sepM <- lapply(seportM, function(x) x %*% tm %*% t(x))

		if(mau == TRUE){
			#Mauchlyの球面性検定
			sep.Lambda <- sapply(sepM, function(x) det(x) / (sum(diag(x)) / nrow(x))^nrow(x))
			sep.m <- 1 - (2 * matdivider^2 + matdivider + 2) / (6 * (cellN - ol) * matdivider)
			sepChi <- -(cellN - ol) * sep.m * log(sep.Lambda)
			sep.df <- ((matdivider * (matdivider + 1)) / 2) - 1
			sep.p <- pchisq(sepChi, ifelse(sep.df == 0, NA, sep.df), lower.tail = FALSE)
		}else if(har == TRUE){
			#Harrisの多標本球面性検定
			sepA <- lapply(seportM, function(x) array(apply(covbuffer, 3, function(y) x %*% y %*% t(x)), dim = c(nrow(x), nrow(x), ol)))

			sepTr <- lapply(sepA, function(y) apply(y, 3, function(x) sum(diag(x))))
			sepSq <- lapply(sepA, function(y) apply(y, 3, function(x) sum(diag(x %*% x))))
			sep.Lambda <- sapply(1:ncol(divpoint), function(x) sum((othN - 1) * sepTr[[x]])^2 / sum((othN - 1) * sepSq[[x]]))
			sepChi <- pmax(0, ((cellN - ol) * matdivider / 2) * ((cellN - ol) * matdivider / sep.Lambda - 1))
			sep.df <- ((ol * matdivider * (matdivider + 1)) / 2) - 1
			sep.p <- pchisq(sepChi, ifelse(sep.df == 0, NA, sep.df), lower.tail = FALSE)
		}else{
			#Mendozaの多標本球面性検定
			sepA <- lapply(seportM, function(x) array(apply(ptm, 3, function(y) x %*% y %*% t(x)), dim = c(nrow(x), nrow(x), ol)))

			sepL1 <- log((cellN-ol)^((cellN-ol)*(matdivider)/2) / sapply(matdivider, function(x) prod((othN-1)^((othN-1) * x / 2))))
			sepL2 <- sapply(sepA, function(y) sum(log(apply(y, 3, function(x) det(x))) * (othN-1)/2))
			sepL3 <- log(sapply(sepA, function(y) sum(diag(apply(y, c(1, 2), function(x) sum(x))/nrow(y))))) * ((cellN-ol) * matdivider / 2)
			sepL <- sepL1 + sepL2 - sepL3
			sep.m <- 1 - ((((cellN-ol) * matdivider^2 * (matdivider + 1) * (2 * matdivider + 1) - (2*(cellN-ol) * matdivider^2)) * 
				sum(1/(othN-1)) - 4) / (6 * (cellN-ol) * matdivider * (ol * matdivider * (matdivider + 1) - 2)))
			sep.m[is.nan(sep.m)] <- 0#NaNが出たところには０を代入
			sepChi <- - 2 * sep.m * sepL
			sep.Lambda <- exp(sepL)
			sep.df <- ((ol * matdivider * (matdivider + 1)) / 2) - 1
			sep.p <- pchisq(sepChi, ifelse(sep.df == 0, NA, sep.df), lower.tail = FALSE)
		}

		#イプシロンを計算する
		sepLB.ep <- 1 / matdivider
		sepGG.ep <- sapply(sepM, function(x) sum(diag(x))^2 / (nrow(x) * sum(x^2)))
		sepHF.ep <- (cellN * matdivider * sepGG.ep - 2) / (matdivider * (cellN - ol - matdivider * sepGG.ep))

		#行列全体での計算結果に追加する
		eps.Lambda <- append(eps.Lambda, sep.Lambda)
		epsChi <- append(epsChi, sepChi)
		eps.df <- append(eps.df, sep.df)
		eps.p <- append(eps.p, sep.p)
		sig.mark <- sig.sign(eps.p)
		LB.ep <- append(LB.ep, sepLB.ep)
		GG.ep <- append(GG.ep, sepGG.ep)
		HF.ep <- append(HF.ep, sepHF.ep)
	}

	#結果をデータフレームにまとめる
	epsitab <- data.frame("Effect" = effect.name, "Dummy" = eps.Lambda, "approx.Chi" = epsChi, "df" = eps.df, 
		"p" = eps.p, "sig.mark" = sig.mark, "LB" = LB.ep, "GG" = GG.ep, "HF" = HF.ep)
	names(epsitab)[2] <- lamlab#ラベルを検定方法に応じたものに変更

	#オプション；IGAのための統計量を計算
	if(iga == TRUE | ciga == TRUE){
		#HuynhのImproved General Approximate Test
		proSj <- lapply(seportM, function(y) array(apply(covbuffer, 3, function(x) y %*% x %*% t(y)), dim = c(nrow(y), nrow(y), ol)))
		proSh <- lapply(seportM, function(y) y %*% (apply(covbuffer, c(1, 2), function(x) sum((1/othN) * x)) / sum(1/othN)) %*% t(y))

		trDt <- sapply(proSh, function(x) sum(diag(x)))
		trDj <- lapply(proSj, function(y) apply(y, 3, function(x) sum(diag(x))))
		trDj2 <- lapply(proSj, function(y) apply(y, 3, function(x) sum(diag(x %*% x))))
		iga.D <- lapply(seportM, function(x) t(x) %*% x)

		iga.h0b <- trDt^2 / sapply(proSh, function(x) sum(diag(x %*% x)))
		iga.b <- ((cellN - ol) * trDt) / sapply(trDj, function(x) sum((othN - 1) * x))

		iga.cue <- array(sapply(1:ol, function(x) ifelse(1:(ol^2) == ol * (x-1) + x, 1, 0)), dim = c(ol, ol, ol))
		iga.Sigstr <- array(rowSums(sapply(1:ol, function(x) iga.cue[,,x] %x% (covbuffer[,,x]/othN[x]))), dim = c(ol*rl, ol*rl))
		iga.g1 <- lapply(iga.D, function(y) array(rowSums(sapply(1:ol, function(x) iga.cue[,,x] %x% (othN[x] * (1-othN[x]/cellN) * y))), dim = c(ol*rl, ol*rl)))
		if(ol == 1){
			iga.h0c <- rep(1, length(matdivider))
			iga.c <- rep(1, length(matdivider))
		}else{
			iga.g2 <- lapply(iga.D, function(z) apply(combn(ol, 2, function(x) array(ifelse(1:(ol^2) == prod(x) | 
				1:(ol^2) == prod(x) + (ol - 1) * abs(diff(x)), 1, 0), dim = c(ol, ol)) %x% (-prod(othN[x]) * z / cellN)), c(1, 2), function(y) sum(y)))
			iga.G <- mapply(function(x, y) x + y, iga.g1, iga.g2, SIMPLIFY = FALSE)
			iga.GS <- lapply(iga.G, function(x) x %*% iga.Sigstr)
			iga.h0c <- sapply(iga.GS, function(x) sum(diag(x))^2) / sapply(iga.GS, function(x) sum(diag(x %*% x)))
			iga.c <- sapply(iga.GS, function(x) ((cellN - ol) * sum(diag(x)))) / sapply(trDj, function(x) ((ol - 1) * sum((othN - 1) * x)))
		}

		iga.h1 <- ifelse(iga.h0b == 1, 1, (cellN * iga.h0b - 2) / (cellN - ol - iga.h0b))
		iga.h2 <- ifelse(iga.h0c == 1, 1, ((ol - 1) * (cellN * iga.h0c - 2 * (ol - 1))) / ((cellN - ol) * (ol - 1) - iga.h0c))
		if(ol == 1){
			iga.eta <- mapply(function(y, z) sum((othN - 1)^3 / ((othN + 1) * (othN - 2)) * (othN * y^2 - 2 * z)), trDj, trDj2)
		}else{
			iga.eta <- mapply(function(y, z) sum((othN - 1)^3 / ((othN + 1) * (othN - 2)) * (othN * y^2 - 2 * z)) + 2 * sum(combn(ol, 2, function(x) prod((othN[x] - 1) * y[x]))), trDj, trDj2)
		}
		iga.sigma <- mapply(function(x, y) sum((othN - 1)^2 / ((othN + 1) * (othN - 2)) * ((othN - 1) * y - x^2)), trDj, trDj2)
		iga.e <- iga.eta / iga.sigma

		#Algina-LecoutreのCorrected Improved General Approximation Testのための指標
		iga.al1 <- ifelse(iga.h0b == 1, 1, ((cellN - ol + 1) * iga.h0b - 2) / (cellN - ol - iga.h0b))
		iga.al2 <- ifelse(iga.h0c == 1, 1, ((ol - 1) * ((cellN - ol + 1) * iga.h0c - 2 * (ol - 1))) / ((cellN - ol) * (ol - 1) - iga.h0c))

		#結果をデータフレームにまとめる
		if(length(effect.name) == 1) iga.name <- effect.name
		else iga.name <- effect.name[-1]

		if(iga == TRUE){
			epsi.info1 <- sub("Epsilons", "Estimates for IGA", epsi.info1)
			if(repleng == 1) epsitab <- cbind(epsitab[,1:6], "b_hat" = iga.b, "c_hat" = iga.c, "h_d" = iga.h1, "h_dd" = iga.h2, "h" = iga.e)
			else epsitab <- cbind(epsitab[,1:6], "b_hat" = c(NA, iga.b), "c_hat" = c(NA, iga.c), "h_d" = c(NA, iga.h1), "h_dd" = c(NA, iga.h2), "h" = c(NA, iga.e))
		}else{
			epsi.info1 <- sub("Epsilons", "Estimates for CIGA", epsi.info1)
			if(repleng == 1) epsitab <- cbind(epsitab[,1:6], "b_hat" = iga.b, "c_hat" = iga.c, "h_d" = iga.al1, "h_dd" = iga.al2, "h" = iga.e)
			else epsitab <- cbind(epsitab[,1:6], "b_hat" = c(NA, iga.b), "c_hat" = c(NA, iga.c), "h_d" = c(NA, iga.al1), "h_dd" = c(NA, iga.al2), "h" = c(NA, iga.e))
		}
	}

	return(list("epsi.info1" = epsi.info1, "epsitab" = epsitab))

}


#分散分析表を作る関数
anova.modeler <- function(dat, design, type2 = FALSE, lb = FALSE, gg = FALSE, hf = FALSE, auto = FALSE, 
	mau = FALSE, har = FALSE, iga = FALSE, ciga = FALSE, peta = FALSE, pomega = FALSE, prep = FALSE, inter = NA){
	#要因計画の型に合わせてss.calc関数を適用し分散分析表を作る
	bet.with <- strsplit(design, "s")
	full.elem <- unlist(lapply(1:nchar(design), function(y) combn(c("s", LETTERS[1:(nchar(design)-1)]), y, function(x) gsub(", ", ":", toString(x)))))

	#効果の自由度
	cellN <- length(unique(dat$s))
	flev <- sapply(2:nchar(design), function(x) length(unique(dat[,x])))#各要因の水準数
	eff.df <- unlist(sapply(1:(nchar(design)-1), function(y) combn(nchar(design)-1, y, function(x) prod(flev[x]-1))))#効果の自由度

	#要因計画のタイプ別の処理
	if(nchar(bet.with[[1]][1]) != 0 && is.na(bet.with[[1]][2])){#被験者間計画の場合
		#full.elemの整理
		er.cue <- gsub(", ", ":", toString(c("s", strsplit(strsplit(design, "s")[[1]][1], "")[[1]])))
		full.elem <- full.elem[-setdiff(grep("s", full.elem), grep(er.cue, full.elem))]
		full.elem <- full.elem[-length(full.elem)]

		#epsilon.calcの出力に対応する情報
		epsi.info1 <- NA
		epsitab <- NA
		ano.info1 <- NA

		#誤差項の自由度
		er.df <- cellN - prod(flev[1:nchar(bet.with[[1]][1])])

	}else if(nchar(bet.with[[1]][1]) != 0){#混合要因計画の場合
		#full.elemの整理
		er.cue <- gsub(", ", ":", toString(c("s", strsplit(strsplit(design, "s")[[1]][1], "")[[1]])))
		full.elem <- full.elem[-setdiff(grep("s", full.elem), grep(er.cue, full.elem))]

		#epsilon.calcの適用
		epsiresults <- epsilon.calc(dat = dat, design = design, mau = mau, har = har, iga = iga, ciga = ciga)
		epsi.info1 <- epsiresults$epsi.info1
		epsitab <- epsiresults$epsitab

		#誤差項の自由度
		er.pri <- cellN - prod(flev[1:nchar(bet.with[[1]][1])])
		er.df <- er.pri * c(1, unlist(sapply(1:nchar(bet.with[[1]][2]), function(y) 
			combn(nchar(bet.with[[1]][2]), y, function(x) prod(flev[-(1:nchar(bet.with[[1]][1]))][x]-1)))))

	}else{#被験者内計画の場合
		#epsilon.calcの適用
		epsiresults <- epsilon.calc(dat = dat, design = design, mau = mau, har = har, iga = iga, ciga = ciga)
		epsi.info1 <- epsiresults$epsi.info1
		epsitab <- epsiresults$epsitab

		#誤差項の自由度
		er.df <- (cellN - 1) * c(1, eff.df)
	}

	#ss.calcの適用
	ss.results <- ss.calc(full.elem = full.elem, dat = dat, type2 = type2)#ss.calc関数にモデルの要素を投入し，平方和を得る

	#並べ替えのための指標と自由度の計算
	eff.elem <- full.elem[-grep("s", full.elem)]#効果の項のみ
	er.elem <- grep("s", full.elem, value = TRUE)#誤差項のみ
	if(length(er.elem) == 0){#被験者間計画の場合
		source.cue <- c(full.elem, "ss.Er", "ss.T")
		df.col <- c(eff.df, er.df, length(dat$y)-1)
		mse.row <- length(source.cue) - 1
	}else{#被験者内要因を含む計画の場合
		source.ord <- sapply(eff.elem, function(x) min(elematch(x, er.elem)))
		source.cue <- unlist(sapply(1:length(er.elem), function(x) c(eff.elem[source.ord == x], er.elem[x])))
		source.cue <- c(source.cue, "ss.T")

		df.col <- unlist(sapply(1:length(er.elem), function(x) c(eff.df[source.ord == x], er.df[x])))
		df.col <- c(df.col, length(dat$y)-1)
		mse.row <- match(er.elem, source.cue)

		#各効果の自由度を調整するための係数をベクトル化する
		if(iga == TRUE){
			mdf <- 1
			ano.info1 <- "== Huynh's Improved General Approximation Test =="
		}else if(ciga == TRUE){
			mdf <- 1
			ano.info1 <- "== Algina-Lecoutre's Corrected Improved General Approximation Test =="
		}else if(lb == TRUE){
			mdf <- pmin(1, rep(c(1, epsitab$LB[epsitab$Effect != "Global"]), c(mse.row[1], diff(mse.row))))
			ano.info1 <- "== Geisser-Greenhouse's Conservative Test =="
		}else if(gg == TRUE){
			mdf <- pmin(1, rep(c(1, epsitab$GG[epsitab$Effect != "Global"]), c(mse.row[1], diff(mse.row))))
			ano.info1 <- "== Adjusted by Greenhouse-Geisser's Epsilon =="
		}else if(hf == TRUE){
			mdf <- pmin(1, rep(c(1, epsitab$HF[epsitab$Effect != "Global"]), c(mse.row[1], diff(mse.row))))
			ano.info1 <- "== Adjusted by Huynh-Feldt's Epsilon =="
		}else if(auto == TRUE){
			sigepsi <- epsitab
			sigepsi$GG[((sigepsi$sig.mark == "") | (sigepsi$sig.mark == "ns"))] <- 1
			mdf <- pmin(1, rep(c(1, sigepsi$GG[sigepsi$Effect != "Global"]), c(mse.row[1], diff(mse.row))))
			ano.info1 <- "== Adjusted by Greenhouse-Geisser's Epsilon for Suggested Violation =="
		}else{
			mdf <- 1
			ano.info1 <- NA
		}

		#自由度を調整する
		df.col <- c(mdf, 1) * df.col
	}

	f.denomi <- rep(mse.row, c(mse.row[1], diff(mse.row)))
	f.denomi[mse.row] <- NA
	f.denomi <- append(f.denomi, NA)

	#分散分析表を作る
	source.col <- gsub("ss.Er", "Error", gsub("ss.T", "Total", gsub(":", "x", source.cue)))
	ss.col <- unlist(sapply(source.cue, function(x) ss.results[names(ss.results) == x]))
	attributes(ss.col) <- NULL#もとの変数名を消去
	ms.col <- ss.col / df.col#MSを計算する
	f.col <- ms.col[1:length(ms.col)] / ms.col[f.denomi]#F値を計算する

	#IGA，CIGAの適用
	if((iga == TRUE && !is.na(epsi.info1[1])) | (ciga == TRUE && !is.na(epsi.info1[1]))){
		mfv <- c(rep(1, mse.row[1]), rep(epsitab$c_hat[epsitab$Effect != "Global"], c(diff(mse.row))))
		mfv[mse.row[-length(mse.row)] + 1] <- epsitab$b_hat[epsitab$Effect != "Global"]
		mfv <- c(mfv, 1)

		iga.df <- c(df.col[1:mse.row[1]], rep(epsitab$h_dd[epsitab$Effect != "Global"], c(diff(mse.row))))
		iga.df[mse.row[-length(mse.row)] + 1] <- epsitab$h_d[epsitab$Effect != "Global"]
		iga.df[mse.row[-1]] <- epsitab$h[epsitab$Effect != "Global"]
		iga.df <- append(iga.df, length(dat$y)-1)

		df.col <- ifelse(iga.df >= df.col, df.col, iga.df)#推定した自由度がもとの自由度よりも高くなった場合は調整を行わない
	}else{
		mfv <- 1#適用しないときは１（調整なし）
	}

	f.col <- f.col / mfv#IGA，CIGAのための推定値によってF値を調整
	p.col <- pf(f.col, df.col, df.col[f.denomi], lower.tail = FALSE)#p値を算出する
	sig.col <- sig.sign(p.col)#p値が有意かどうかを判定して記号を表示する
	anovatab <- data.frame(source.col, ss.col, df.col, ms.col, f.col, p.col, sig.col)#分散分析表をまとめたデータフレーム

	#効果量の計算と追加
	if(peta == TRUE){#偏イータ二乗
		peta.col <- ss.col / (ss.col + ss.col[f.denomi])
		anovatab <- cbind(anovatab, "p.eta^2" = peta.col)
	}
	if(pomega == TRUE){#偏オメガ二乗；値が負になったときは０に直す
		rscue <- c(mse.row, length(source.col))
		eachN <- rep(NA, length(source.col))
		eachN[-rscue] <- sapply(source.col[-rscue], function(y) sum(tapply(dat$s, dat[strsplit(y, "x")[[1]]], function(x) length(unique(x)))))
		pomega.col <- pmax(0, (ss.col - df.col * ms.col[f.denomi]) / (ss.col + (eachN - df.col) * ms.col[f.denomi]))
		anovatab <- cbind(anovatab, "p.omega^2" = pomega.col)
	}
	if(prep == TRUE){#p_rep（両側）
		prep.col <- pmin(0.9999, pnorm(qnorm(1 - p.col/2) / sqrt(2)))
		anovatab <- cbind(anovatab, "p_rep" = prep.col)
	}

	#結果を返す；interに入力がない場合はanovatabとmse.row，入力がある場合はintertabを返す
	if(is.na(inter)){
		return(list("epsi.info1" = epsi.info1, "epsitab" = epsitab, "ano.info1" = ano.info1, "anovatab" = anovatab, "mse.row" = mse.row))
	}else{
		#単純主効果の検定用の出力を用意する
		sim.row <- charmatch(inter, anovatab$source.col)
		intertab <- rbind(anovatab[sim.row,], anovatab[f.denomi[sim.row],])

		#球面性検定の結果からinterに関する部分のみ取り出す
		if(charmatch(inter, strsplit(design, "")[[1]]) < charmatch("s", strsplit(design, "")[[1]])){
			#被験者間要因ならNAの行を返す
			if(iga == TRUE | ciga == TRUE) interepsi <- rep(NA, 11)#IGA，CIGAを使ったときは列数が多い
			else interepsi <- rep(NA, 9)
		}else{
			interepsi <- epsitab[charmatch(inter, epsitab$Effect),]
		}
		return(list("intertab" = intertab, "interepsi" = interepsi))
	}
}


#修正Bonferroniの方法（Holmの方法，Shafferの方法，Holland-Copenhaverの方法）による多重比較を行う関数
#デフォルトはShafferの方法を適用し，holm = TとするとHolmの方法，hc = TとするとHolland-Copenhaverの方法を適用する
#s2r = T，s2d = Tとすると，具体的な棄却のパターンを反映した有意水準の調整によるShafferの方法を適用する
mod.Bon <- function(dat, design, taref, bet.mse, factlabel = NA, type2 = FALSE, holm = FALSE, hc = FALSE, 
	s2r = FALSE, s2d = FALSE, fs1 = FALSE, fs2r = FALSE, fs2d = FALSE, alpha = 0.05, criteria = FALSE){
	#対象となる要因のラベルを得る
	if(is.na(factlabel)) factlabel <- taref

	bonl <- nlevels(dat[,taref])#水準数の取得
	h0size <- bonl * (bonl-1)/2#帰無仮説の個数
	bon.name <- combn(bonl, 2, function(ij) paste(levels(dat[,taref])[ij][1], "-", levels(dat[,taref])[ij][2], sep = ""))#可能な対の組み合わせのラベル
	bon.num <- charmatch(taref, names(dat))-1#分析対象となる要因が何番目の要因かを特定

	#周辺平均を計算する
	cont.means <- tapply(dat$y, dat[,2:nchar(design)], mean)#各セルの平均を求める
	bon.means <- apply(cont.means, bon.num, mean)#分析対象となる周辺平均を求める

	factlevels <- sapply(names(dat), function(x) nlevels(dat[,x]))#各要因の水準数
	factlevels[charmatch(taref, names(dat))] <- 2#多重比較の対象となる効果の水準数を２に固定
	factlevels <- factlevels[!(factlevels == factlevels[1] | factlevels == factlevels[length(factlevels)])]#最初と最後（sとyの列）を除く

	cont.N <- table(dat[,2:nchar(design)])#各セルのデータ数
	bon.denomi <- apply(1/cont.N, bon.num, mean) / (prod(factlevels)/2)#セルごとに重み付けしたデータ数の平均を残りの条件数で割ったもの

	bon.delta <- combn(bonl, 2, function(ij) (bon.means[ij][1] - bon.means[ij][2]))#平均偏差

	#分析対象が被験者間要因か被験者内要因かによって標準誤差を得る方法を変える
	if(length(bet.mse) != 1){
		#被験者間要因の場合；上位の分析のMSeを適用する
		bon.df <- bet.mse$df.col#自由度
		bon.Ve <- bet.mse$ms.col#平均平方
	}else{
		#被験者内要因の場合；比較する２水準ごとにMSeを再計算する
		bon.lev <- combn(levels(dat[,taref]), 2, function(x) x)#各水準の組み合わせ
		subdat <- apply(bon.lev, 2, function(x) dat[dat[, taref] == x[1] | dat[, taref] == x[2],])#データを多重比較の対象となる効果について２水準ずつのサブデータにリスト化
		for(i in 1:length(subdat)){
			subdat[[i]][, taref] <- subdat[[i]][, taref][, drop = TRUE]
		}

		bon.anova <- lapply(subdat, function(x) anova.modeler(dat = x, design = design, type2 = type2, 
			inter = taref)$intertab[2,])
		bon.df <- sapply(bon.anova, function(x) x$df.col)#自由度
		bon.Ve <- sapply(bon.anova, function(x) x$ms.col)
	}

	#検定統計量とｐ値を得る
	bon.SE <- sqrt(combn(bonl, 2, function(ij) sum(bon.denomi[ij])) * bon.Ve)
	bon.t <- abs(bon.delta / bon.SE)#ｔ値は絶対値を取る
	bon.p <- pt(bon.t, bon.df, lower.tail = FALSE) * 2#両側確率

	#結果をデータフレームにまとめる
	bontab <- data.frame("pair" = bon.name, "interval" = bon.delta, "t" = bon.t, "df" = bon.df, "p.value" = bon.p)
	bontab <- bontab[order(bontab$p.value),]#p値の小さい順に並べ替え

	#調整した有意水準を設定する
	if(holm == TRUE){
		p.criteria <- h0size:1#Holmの方法用の調整値
		bon.info1 <- paste("== Holm's Sequentially Rejective Bonferroni Procedure ==", sep = "")

	}else if(s2d | fs2d == TRUE){#Donoguhe（2004）のアルゴリズムをベースとするShafferの多重比較のための論理ステップの計算
		#Donoghue, J. R. (2004). Implementing Shaffer's multiple comparison procedure for a large number of groups.
		#Recent developments in multiple comparison procedures (Institute of mathematical statistics-Lecture Notes-Monograph Series, 47), pp. 1-23.
		#Donoghueと完全に同じ手順ではないことに注意

		#隣接行列を作る
		bon.comb <- combn(bonl, 2, function(x) x)#帰無仮説を表す行列
		bon.comb <- bon.comb[,order(bon.p)]#ｐ値の順に並べ替え
		hvec <- 1:bonl
		a.mat <- diag(bonl)
		diag(a.mat) <- 0

		shaf.value <- c(h0size, rep(NA, h0size - 1))#ステップごとの仮説数を代入するためのベクトル

		allcomb <- lapply((bonl-1):1, function(y) combn(bonl, y, function(x) x, simplify = FALSE))#すべての帰無仮説の組み合わせを示すリスト
		allcomb <- unlist(allcomb, recursive = FALSE)#リストの階層をなくす

		for(j in 1:(h0size-1)){
			#隣接行列に棄却された仮説を書き込む
			a.mat[bon.comb[1, j], bon.comb[2, j]] <- 1#棄却された帰無仮説の部分に１を代入
			a.mat[bon.comb[2, j], bon.comb[1, j]] <- 1#対角線を通して反対の側にも代入

			#未分化クラスを作る
			#隣接行列の下位行列の中から０のみで構成される正方行列を探す
			#未分化クラスを表す行列：各行が各水準に相当；各列が帰無仮説を表す（互いに差がない水準に１を代入）
			undiff <- array(rep(0, bonl), c(bonl, 1))#ダミー
			cnt <- 1
			while(cnt <= length(allcomb)){
				hnum <- allcomb[[cnt]]
				if(max(colSums(undiff[hnum, , drop = FALSE])) == length(hnum)){
				#上位の仮説に包含される仮説は含めない；カットはしない
					cnt <- cnt + 1
				}else if(sum(a.mat[hnum, hnum]) == 0){
				#正方行列の場合は成立する帰無仮説を表す列を追加
					undiff <- cbind(undiff, 1 - 0^match(hvec, hnum, nomatch = 0))
					cnt <- cnt + 1
				}else{
				#その他の場合；このパターンは後に支持されることはないので，allcombからカット
					allcomb <- allcomb[-cnt]
					}
				}
			undiff <- undiff[,-1]#一列目のダミーを除く
			gsize <- colSums(undiff)#各グループの要素数を示すベクトル

			#sig.minを決定する
			sig.min <- max(gsize)^2#最大クラスの要素数を二乗した値
			nxcand <- undiff#未分化クラスのコピー
			gi <- 1
			while(ncol(nxcand) > 1 && nrow(nxcand) > 1){
				nxcand <- nxcand[(1:nrow(nxcand))[nxcand[,gi] == 0], , drop = FALSE]
				lengvec <- colSums(nxcand)#各クラスの要素数
				sig.min <- sig.min + max(lengvec)^2#最大クラスの要素数の二乗値を足す
				gi <- which.max(lengvec)#最大クラスの番号
				}

			#don.maxを決定する
			don.smax <- sig.min

			for(i in 2:min(ncol(undiff)-1, bonl)){
				don.sig <- gsize[i]^2
				nxcand <- undiff
				gi <- i
				while(ncol(nxcand) > 1 && nrow(nxcand) > 1){
					nxcand <- nxcand[(1:nrow(nxcand))[nxcand[,gi] == 0], , drop = FALSE]
					lengvec <- colSums(nxcand)#各クラスの要素数
					don.sig <- don.sig + max(lengvec)^2#最大クラスの要素数の二乗値を足す
					gi <- which.max(lengvec)#最大クラスの番号
					}
				don.smax <- max(don.sig, don.smax)#より大きい値を残す
				}

			shaf.value[j+1] <- (don.smax - bonl) / 2

			}

		if(s2d == T){
			p.criteria <- shaf.value#Shafferの方法用の調整値
			bon.info1 <- paste("== Shaffer's Modified Sequentially Rejective Bonferroni Procedure [SPECIFIC] ==", "\n", 
				"== This computation is based on the algorithm by Donoghue (2004). ==", sep = "")
			shaf.meth <- paste(" [SPECIFIC] ==", "\n", "== This computation is based on the algorithm by Donoghue (2004). ==", sep = "")
		}else{
			p.criteria <- c(shaf.value[2], shaf.value[2:length(shaf.value)])#F-Shafferの方法用の調整値
			bon.info1 <- paste("== Shaffer's F-Modified Sequentially Rejective Bonferroni Procedure [SPECIFIC] ==", "\n", 
				"== This computation is based on the algorithm by Donoghue (2004). ==", sep = "")
			shaf.meth <- paste(" [SPECIFIC] ==", "\n", "== This computation is based on the algorithm by Donoghue (2004). ==", sep = "")
		}

	}else{#Rasmussen（1993）のアルゴリズムによるShafferの多重比較のための論理ステップの計算
		#Rasmussen, J. L. (1993). Algorithm for Shaffer's multiple comparison tests. Educational and Psychological Measurement, 53, 329-335.

		#平均間の異同パターンを表す行列を作る
		hpattern <- 2^(bonl-1)#可能な真偽の仮説のパターン数
		nbuffer <- 2^((bonl-2):0)
		c.mat <- cbind(rep(0, hpattern), sapply(nbuffer, function(x) rep(c(rep(0, x), rep(1, x)), nbuffer[1]/x)))
		c.mat <- t(apply(c.mat, 1, function(x) cumsum(x)))

		f.mat <- combn(bonl, 2, function(x) c.mat[,x[1]] - c.mat[,x[2]])#各水準の組み合わせを表現する行列
		f.mat[f.mat != 0] <- 1#帰無仮説が真のときに０，偽のときに１となるようにする
		rebon.p <- bon.p[order(combn(rank(bon.means), 2, function(x) prod(x)))]#ｐ値の順序を平均値の大きさにそって並べ替え
		f.mat <- f.mat[,order(rebon.p)]#ｐ値の小さい順に列を並べ替え
		i.vector <- rowSums(f.mat)#棄却される帰無仮説の数
		t.vector <- h0size - i.vector#成立しうる真の帰無仮説の数

		if(s2r | fs2r == TRUE){#各比較までの特定の仮説が偽であったときの可能な真の帰無仮説の最大数
			shaf.value <- c(max(t.vector), max(t.vector[i.vector >= (2 - 1)][(f.mat[i.vector >= (2 - 1), 1:(2-1)]) == (2 - 1)]))
			shaf.value <- c(shaf.value, sapply(3:h0size, function(x) max(t.vector[i.vector >= (x - 1)][rowSums(f.mat[i.vector >= (x - 1), 1:(x-1)]) == (x - 1)])))
			shaf.meth <- paste(" [SPECIFIC] ==", "\n", "== This computation is based on the algorithm by Rasmussen (1993). ==", sep = "")
		}else{#各比較までの任意の仮説が偽であったときの可能な真の帰無仮説の最大数
			shaf.value <- sapply(1:h0size, function(x) max(t.vector[i.vector >= (x - 1)]))
			shaf.meth <- " =="
		}

		if(fs1 | fs2r == TRUE){
			p.criteria <- c(shaf.value[2], shaf.value[2:length(shaf.value)])#F-Shafferの方法用の調整値
			bon.info1 <- paste("== Shaffer's F-Modified Sequentially Rejective Bonferroni Procedure", shaf.meth, sep = "")
		}else{
			p.criteria <- shaf.value#Shafferの方法用の調整値
			bon.info1 <- paste("== Shaffer's Modified Sequentially Rejective Bonferroni Procedure", shaf.meth, sep = "")
		}
	}

	#平均値の差の方向を調べ，不等号のベクトルを作る
	bon.differ <- ifelse(bontab$interval <= 0, sub("-", " < ", bontab$pair), sub("-", " > ", bontab$pair))
	#差が見られなかった場合の等号のベクトルを作る
	bon.equal <- sub("-", " = ", bontab$"pair")

	if(criteria == TRUE){#データフレームに調整済み有意水準の列を加える
		if(hc == TRUE){bontab <- transform(bontab, "criteria" = 1 - (1 - alpha) ^ (1/p.criteria))#Sidakの不等式による有意水準の調整
			if(holm == TRUE) bon.info1 <- paste("== Holm's Sequentially Rejective Sidak Procedure ==", sep = "")
			else bon.info1 <- paste("== Holland-Copenhaver's Improved Sequentially Rejective Sidak Procedure", shaf.meth, sep = "")

			if(length(bet.mse) == 1) bon.info1 <- append(bon.info1, "*** CAUTION! This procedure might be inappropriate for dependent means. ***")
		}else{bontab <- transform(bontab, "criteria" = alpha/p.criteria)#Bonferroniの不等式による有意水準の調整
		}
		#有意であった行は不等号，そうでない行は等号を表示する
		bon.sign <- ifelse(cummin(bontab$p.value < bontab$criteria), paste(bon.differ, "*", sep = " "), paste(bon.equal, " ", sep = " "))
	}else{#データフレームに調整済みｐ値の列を加える
		if(hc == TRUE){bontab <- transform(bontab, "adj.p" = pmin(1, cummax((1-(1-bontab$p.value)^p.criteria))))#Sidakの不等式による調整済みｐ値
			if(holm == TRUE) bon.info1 <- paste("== Holm's Sequentially Rejective Sidak Procedure ==", sep = "")
			else bon.info1 <- paste("== Holland-Copenhaver's Improved Sequentially Rejective Sidak Procedure", shaf.meth, sep = "")

			if(length(bet.mse) == 1) bon.info1 <- append(bon.info1, "*** CAUTION! This procedure might be inappropriate for dependent means. ***")
		}else{bontab <- transform(bontab, "adj.p" = pmin(1, cummax(bontab$p.value * p.criteria)))#Bonferroniの不等式による調整済みｐ値
		}
		#有意であった行は不等号，そうでない行は等号を表示する
		bon.sign <- ifelse(bontab$adj.p < alpha, paste(bon.differ, "*", sep = " "), paste(bon.equal, " ", sep = " "))
	}

	#判定結果をデータフレームに反映する
	bontab <- transform(bontab, "significance" = bon.sign)

	#記述統計量の計算
	b.sncol <- tapply(dat$y, dat[,taref], length)#セルごとのデータ数を計算
	b.sdcol <- tapply(dat$y, dat[,taref], sd)#セルごとの標準偏差を計算

	bonstat <- data.frame("Dummy" = levels(dat[, taref]), "N" = b.sncol, "Mean" = bon.means, "S.D." = b.sdcol)

	names(bonstat)[1] <- factlabel#水準を表すラベルをfactlabelとして入力した値に置き換える

	#その他の出力を準備する
	bon.info2 <- if(length(bet.mse) != 1) paste("== The factor < ", factlabel, " > is analysed as independent means. ==", sep  = "")
		else paste("== The factor < ", factlabel, " > is analysed as dependent means. ==", sep  = "")
	bon.info3 <- paste("== Alpha level is ", alpha, ". ==", sep  = "")

	return(list(factlabel, bon.info1, bon.info2, bon.info3, bonstat, bontab))

}


#下位検定を行う関数
post.analyses <- function(dat, design, mainresults, type2 = FALSE, holm = FALSE, hc = FALSE, 
	s2r = FALSE, s2d = FALSE, fs1 = FALSE, fs2r = FALSE, fs2d = FALSE, criteria = FALSE, 
	lb = FALSE, gg = FALSE, hf = FALSE, auto = FALSE, mau = FALSE, har = FALSE, iga = FALSE, ciga = FALSE, 
	peta = FALSE, pomega = FALSE, prep = FALSE){
	anovatab <- mainresults$anovatab

	#要因計画の型から被験者間要因と被験者内要因の情報を得る
	bet.with <- strsplit(design, "s")
	#効果が有意であった行のソースラベルを得る
	sig.source <- anovatab$"source.col"[!((anovatab$"sig.col" == "") | (anovatab$"sig.col" == "ns"))]

	if(length(sig.source) == 0) {
		return(NA)#有意な行がなければここで終了
	}else{
		#下位検定の結果を格納するための空のデータフレームを宣言
		postresults <- lapply(paste("post", 1:length(sig.source), sep = ""), function(x) x)

		#pro.fraction関数を反復適用
		for (i in 1:length(sig.source)) postresults[[i]] <- pro.fraction(dat = dat, design = design, postplan = sig.source[i], 
			bet.with = bet.with, mainresults = mainresults, type2 = type2, holm = holm, hc = hc, 
			s2r = s2r, s2d = s2d, fs1 = fs1, fs2r = fs2r, fs2d = fs2d, criteria = criteria, 
			lb = lb, gg = gg, hf = hf, auto = auto, mau = mau, har = har, iga = iga, ciga = ciga, peta = peta, pomega = pomega, prep = prep)

		return(postresults)
	}
}


#効果のタイプに適した下位検定を割り当てる関数
pro.fraction <- function(dat, design, postplan, bet.with, mainresults, type2 = FALSE, holm = FALSE, hc = FALSE, 
	s2r = FALSE, s2d = FALSE, fs1 = FALSE, fs2r = FALSE, fs2d = FALSE, criteria = FALSE, 
	lb = FALSE, gg = FALSE, hf = FALSE, auto = FALSE, mau = FALSE, har = FALSE, iga = FALSE, ciga = FALSE, 
	peta = FALSE, pomega = FALSE, prep = FALSE){
	#情報の展開
	anovatab <- mainresults$anovatab
	mse.row <- mainresults$mse.row

	#ソースラベルを文字列に変換し，その文字数を得る
	sig.term <- as.character(postplan)
	sig.num <- nchar(sig.term)

	#効果の種類によって違った処理を割り当てる
	if(sig.num > 3){
	#高次の交互作用：効果のラベルを返す
		return(sig.term)

	}else if(sig.num == 3){
	#１次の交互作用については，単純主効果の検討を行う
		each.term <- strsplit(sig.term, "x")#交互作用の各項を分離する
		interfact <- as.numeric(sapply(each.term[[1]], function(x) grep(x, bet.with[[1]][1])))#被験者間要因なら１，被験者内要因なら０を返す

		#交互作用を構成する第一項の項の列番号，水準数を取得
		col.num1 <- charmatch(each.term[[1]][1], names(dat))#列番号
		level.num1 <- nlevels(dat[, col.num1])#水準数

		#交互作用を構成する第二項の項の列番号，水準数を取得
		col.num2 <- charmatch(each.term[[1]][2], names(dat))#列番号
		level.num2 <- nlevels(dat[, col.num2])#水準数

		#単純主効果の検定用の要因計画を作る
		resdesign1 <- sub(each.term[[1]][2], "", design)#つぶす要因のラベルを消去
		bet.with1 <- strsplit(resdesign1, "s")#残りの要因を被験者間と被験者内に分離
		bet.num1 <- if(nchar(bet.with1[[1]][1]) == 0) 0 else 1:nchar(bet.with1[[1]][1])#被験者間要因の数を数える
		with.num1 <- if(is.na(bet.with1[[1]][2])) 0 else (max(bet.num1)+1):(nchar(resdesign1)-1)#被験者内要因の数を数える
		bet1 <- gsub(", ", "", toString(LETTERS[bet.num1]))#被験者間要因のラベルを新たに作成
		with1 <- gsub(", ", "", toString(LETTERS[with.num1]))#被験者内要因のラベルを新たに作成
		postdesign1 <- paste(bet1, "s", with1, sep = "")#両要因を合成

		resdesign2 <- sub(each.term[[1]][1], "", design)#つぶす要因のラベルを消去
		bet.with2 <- strsplit(resdesign2, "s")#残りの要因を被験者間と被験者内に分離
		bet.num2 <- if(nchar(bet.with2[[1]][1]) == 0) 0 else 1:nchar(bet.with2[[1]][1])#被験者間要因の数を数える
		with.num2 <- if(is.na(bet.with2[[1]][2])) 0 else (max(bet.num2)+1):(nchar(resdesign2)-1)#被験者内要因の数を数える
		bet2 <- gsub(", ", "", toString(LETTERS[bet.num2]))#被験者間要因のラベルを新たに作成
		with2 <- gsub(", ", "", toString(LETTERS[with.num2]))#被験者内要因のラベルを新たに作成
		postdesign2 <- paste(bet2, "s", with2, sep = "")#両要因を合成

		#分析対象となる効果が被験者間要因か被験者内要因かによってソース列のラベルを作成
		if(is.na(interfact[1]) == FALSE) err.label1 <- "Er"
		else err.label1 <- paste("sx", each.term[[1]][1], sep = "")

		if(is.na(interfact[2]) == FALSE) err.label2 <- "Er"
		else err.label2 <- paste("sx", each.term[[1]][2], sep = "")

		#データフレームを各要因の各水準ごとに分離
		subdat1 <- dat#データフレームをコピー
		rename.vector1 <- c("s", LETTERS[1:(nchar(design)-2)], "y")
		rename.vector1 <- append(rename.vector1, "X", after = col.num2-1)
		names(subdat1) <- rename.vector1#列名を変更
		target.effect1 <- names(subdat1)[col.num1]#検討したい効果の変更後の列名を取得
		subdat1 <- split(subdat1, subdat1[,col.num2])#分析対象でない要因の水準ごとにデータフレームを分割
		subdat1 <- lapply(subdat1, function(x) x[-col.num2])#X列を消去

		subdat2 <- dat#データフレームをコピー
		rename.vector2 <- c("s", LETTERS[1:(nchar(design)-2)], "y")
		rename.vector2 <- append(rename.vector2, "X", after = col.num1-1)
		names(subdat2) <- rename.vector2#列名を変更
		target.effect2 <- names(subdat2)[col.num2]#検討したい効果の変更後の列名を取得
		subdat2 <- split(subdat2, subdat2[,col.num1])#分析対象でない要因の水準ごとにデータフレームを分割
		subdat2 <- lapply(subdat2, function(x) x[-col.num1])#X列を消去

		#各効果について要因を１つ落とした分散分析
		sim.effects1 <- rep(list(NA), col.num2)#結果格納用リスト
		sim.effects2 <- rep(list(NA), col.num1)#結果格納用リスト
		sim.effects1 <- lapply(subdat1, function(x) anova.modeler(dat = x, design = postdesign1, type2 = type2, 
			lb = lb, gg = gg, hf = hf, auto = auto, mau = mau, har = har, iga = iga, ciga = ciga, 
			peta = peta, pomega = pomega, prep = prep, inter = target.effect1))
		sim.effects2 <- lapply(subdat2, function(x) anova.modeler(dat = x, design = postdesign2, type2 = type2, 
			lb = lb, gg = gg, hf = hf, auto = auto, mau = mau, har = har, iga = iga, ciga = ciga, 
			peta = peta, pomega = pomega, prep = prep, inter = target.effect2))

		#結果を１つのデータフレームにまとめる
		simtab <- data.frame()#結果格納用データフレーム
		simepsi <- data.frame()#結果格納用データフレーム
		subsource <- c()#ソースラベル格納用ベクトル

		#計画のタイプによって出力を整理し直す
		if(substr(design, nchar(design), nchar(design)) == "s"){
			#被験者間計画なら，anovatabから主分析の誤差平方和等を得る
			between.mse <- anovatab[charmatch("Error", anovatab$source.col),]
			bet.mse1 <- between.mse
			bet.mse2 <- between.mse

			#主分析の誤差平方和を用いて単純主効果の結果を再計算
			for (i in 1:level.num2){
				simtab <- rbind(simtab, sim.effects1[[i]]$intertab[1,])
				subsource <- append(subsource, c(paste(each.term[[1]][1], " at ", names(sim.effects1)[i], sep = "")))
			}
			for (i in 1:level.num1){
				simtab <- rbind(simtab, sim.effects2[[i]]$intertab[1,])
				subsource <- append(subsource, c(paste(each.term[[1]][2], " at ", names(sim.effects2)[i], sep = "")))
			}
			simtab <- rbind(simtab, between.mse)
			sim.f.col <- simtab$ms.col[1:(nrow(simtab)-1)] / simtab$ms.col[nrow(simtab)]
			sim.p.col <- pf(sim.f.col[1:(nrow(simtab)-1)], simtab$df.col[1:(nrow(simtab)-1)], simtab$df.col[nrow(simtab)], lower.tail = FALSE)
			sim.sig.col <- sig.sign(sim.p.col)

			#計算結果をデータフレームに反映
			simtab$sig.col <- ""
			simtab <- cbind(simtab[,1:4], "f.col" = c(sim.f.col, NA), "p.col" = c(sim.p.col, NA), "sig.col" = c(sim.sig.col, ""))

			#効果量の再計算
			if(peta == TRUE){#偏イータ二乗
				peta.col <- simtab$ss.col[1:(nrow(simtab)-1)] / (simtab$ss.col[1:(nrow(simtab)-1)]  + simtab$ss.col[nrow(simtab)])
				simtab <- cbind(simtab, "p.eta^2" = c(peta.col, NA))
			}
			if(pomega == TRUE){#偏オメガ二乗
				pomega.col <- pmax(0, (simtab$ss.col[1:(nrow(simtab)-1)] - simtab$df.col[1:(nrow(simtab)-1)] * simtab$ms.col[nrow(simtab)]) / 
					(simtab$ss.col[1:(nrow(simtab)-1)] + (cellN - simtab$df.col[1:(nrow(simtab)-1)]) * simtab$ms.col[nrow(simtab)]))
				simtab <- cbind(simtab, "p.omega^2" = c(pomega.col, NA))
			}
			if(prep == TRUE){#p_rep（両側）
				prep.col <- pmin(0.9999, pnorm(qnorm(1 - sim.p.col/2) / sqrt(2)))
				simtab <- cbind(simtab, "p_rep" = c(prep.col, NA))
			}

			subsource <- append(subsource, "Error")
			simepsi <- NA#球面性検定の結果はなし

			#ソース列，行番号のラベルを張り替える
			simtab$source.col <- subsource
			row.names(simtab) <- c(1:nrow(simtab))
		}else{
			#被験者内計画，混合要因計画の場合，結果を１行ずつ取り出してrbindで結合する
			for (i in 1:level.num2){
				simtab <- rbind(simtab, sim.effects1[[i]]$intertab)
				simepsi <- rbind(simepsi, sim.effects1[[i]]$interepsi)
				subsource <- append(subsource, c(paste(each.term[[1]][1], " at ", names(sim.effects1)[i], sep = ""), 
					paste(err.label1, " at ", names(sim.effects1)[i], sep = "")))
			}
			names(simepsi) <- names(sim.effects2[[1]]$interepsi)#名前が異なるとrbindできないので統一；混合要因計画の場合，必ず被験者間要因がsim.effects1に入ることを前提とする
			for (i in 1:level.num1){
				simtab <- rbind(simtab, sim.effects2[[i]]$intertab)
				simepsi <- rbind(simepsi, sim.effects2[[i]]$interepsi)
				subsource <- append(subsource, c(paste(each.term[[1]][2], " at ", names(sim.effects2)[i], sep = ""), 
					paste(err.label2, " at ", names(sim.effects2)[i], sep = "")))
			}
			simepsi$Effect <- subsource[(1:length(subsource)) %% 2 == 1]#subsourceの奇数番の値をsimepsiの効果ラベルに貼り付ける
			simepsi <- simepsi[!is.na(simepsi$df),]#dfがNAの行（被験者間効果の行）を除く
			if(nrow(simepsi) == 0) simepsi <- NA#数値のある行がなければNAを代入

			#ソース列，行番号のラベルを張り替える
			simtab$source.col <- subsource
			row.names(simtab) <- c(1:nrow(simtab))

			#混合計画における被験者間要因に対して主分析の平均平方を取り出す
			if(substr(postdesign1, nchar(postdesign1), nchar(postdesign1)) == "s"){
				bet.mse1 <- simtab[1:level.num2 * 2,]#被験者間要因が分析対象の場合
			}else{
				bet.mse1 <- NA#被験者内要因が分析対象の場合
			}

			if(substr(postdesign2, nchar(postdesign2), nchar(postdesign2)) == "s"){
				bet.mse2 <- simtab[(level.num2 + 1):(level.num1 + level.num2) * 2,]#被験者間要因が分析対象の場合		
			}else{
				bet.mse2 <- NA#被験者内要因が分析対象の場合
			}
		}

		#記述統計量の計算
		sim.sncol <- as.vector(tapply(dat$y, list(dat[,col.num2], dat[,col.num1]), length))#セルごとのデータ数を計算
		sim.mncol <- as.vector(tapply(dat$y, list(dat[,col.num2], dat[,col.num1]), mean))#セルごとの平均を計算
		sim.sdcol <- as.vector(tapply(dat$y, list(dat[,col.num2], dat[,col.num1]), sd))#セルごとの標準偏差を計算

		sim.stat <- data.frame("Term1" = rep(levels(dat[,col.num1]), each = level.num2), 
			"Term2" = rep(levels(dat[,col.num2]), level.num1), "N" = sim.sncol, "Mean" = sim.mncol, "S.D." = sim.sdcol)

		#行ラベルの張り替え
		names(sim.stat)[1] <- each.term[[1]][1]
		names(sim.stat)[2] <- each.term[[1]][2]

		#有意であった行のソースをチェック
		sim.sig.col <- simtab$source.col[!((simtab$sig.col == "") | (simtab$sig.col == "ns"))]

		#多重比較を実行
		if(length(sim.sig.col) == 0){
			sim.multresults <- NA#有意な行がない場合はNAを代入
		}else{
			#多重比較の結果を格納するための空のデータフレームを宣言
			sim.multresults <- sapply(paste("mult", 1:length(sim.sig.col), sep = ""), function(x) list(x))

			#多重比較の関数を反復適用
			for (i in 1:length(sim.sig.col)){
				#ソースラベルを分解
				source.term <- strsplit(sim.sig.col[i], " at ")

				if((source.term[[1]][1] == each.term[[1]][1]) && (level.num1 >= 3) == TRUE){
					#第一項の単純主効果が有意；subdat1を使用
					multdat <- subdat1[[source.term[[1]][2]]]#対象となる水準のデータのみ抽出
					#混合要因計画の被験者間要因の場合は複数の行から適切なMSeを選択；それ以外の場合は単一の行をコピー
					if(!is.na(bet.mse1) && nrow(bet.mse1) != 1) bet.mse <- bet.mse1[grep(source.term[[1]][2], bet.mse1$source.col),]
					else bet.mse <- bet.mse1
					sim.multresults[[i]] <- mod.Bon(dat = multdat, design = postdesign1, taref = target.effect1, bet.mse = bet.mse, factlabel = sim.sig.col[i], 
						type2 = type2, holm = holm, hc = hc, s2r = s2r, s2d = s2d, fs1 = fs1, fs2r = fs2r, fs2d = fs2d, criteria = criteria)
				}else if((source.term[[1]][1] == each.term[[1]][2]) && (level.num2 >= 3) == TRUE){
					#第二項の単純主効果が有意；subdat2を使用
					multdat <- subdat2[[source.term[[1]][2]]]#対象となる水準のデータのみ抽出
					#混合要因計画の被験者間要因の場合は複数の行から適切なMSeを選択；それ以外の場合は単一の行をコピー
					if(!is.na(bet.mse2) && nrow(bet.mse2) != 1) bet.mse <- bet.mse2[grep(source.term[[1]][2], bet.mse2$source.col),]
					else bet.mse <- bet.mse2
					sim.multresults[[i]] <- mod.Bon(dat = multdat, design = postdesign2, taref = target.effect2, bet.mse = bet.mse, factlabel = sim.sig.col[i], 
						type2 = type2, holm = holm, hc = hc, s2r = s2r, s2d = s2d, fs1 = fs1, fs2r = fs2r, fs2d = fs2d, criteria = criteria)
				}else{sim.multresults[[i]] <- NA}
			}
		}

		return(list("sig.term" = sig.term, "sim.stat" = sim.stat, "simepsi" = simepsi, "simtab" = simtab, "sim.multresults" = sim.multresults))

	}else if(sig.num == 1){
	#単純主効果が有意で，水準数が３以上であれば多重比較を行う
		#分析対象となるの項の列番号，水準数を取得
		col.num <- charmatch(sig.term, names(dat))#列番号
		level.num <- nlevels(dat[, col.num])#水準数

		#水準数が３以上の場合にのみ，mod.Bon関数を適用する
		if(level.num >= 3){
			multfact <- grep(sig.term, bet.with[[1]][1])#被験者間要因なら１，被験者内要因なら０を返す

			#被験者間要因か被験者内要因かを判定して，mod.Bon関数にデータを送る
			if(is.na(multfact[1]) == FALSE){
				if(substr(design, nchar(design), nchar(design)) == "s"){
				#被験者間計画なら全体の誤差項のMSを得る
					bet.mse <- anovatab[charmatch("Error", anovatab$source.col),]
				}else{
				#混合要因計画ならその効果の誤差項のMSを得る
					f.denomi <- rep(mse.row, c(mse.row[1], diff(mse.row)))
					bet.mse <- anovatab[f.denomi[charmatch(sig.term, anovatab$source.col)],]
				}
			}else{
				bet.mse <- NA
			}

			bonout <- mod.Bon(dat = dat, design = design, taref = sig.term, bet.mse = bet.mse, factlabel = NA, 
				type2 = type2, holm = holm, hc = hc, s2r = s2r, s2d = s2d, fs1 = fs1, fs2r = fs2r, fs2d = fs2d, criteria = criteria)

			return(bonout)

		}else{
			return(NA)#水準数が２ならNAを返す
		}
	}
}


#基本情報を出力する関数
info.out <- function(info1, info2, info3){
	cat("\n")#１行空ける
	cat(sprintf(info1), "\n", "\n")
	cat(sprintf(info2), "\n")
	cat(sprintf(info3), "\n", "\n", "\n")
}


#平均と標準偏差の表を出力する関数
bstat.out <- function(bstatist, maxfact, margin){
	cat(sprintf("<< DESCRIPTIVE STATISTICS >>"), sep = "", "\n", "\n")#タイトル
	cat(sprintf(rep("-", (36 + 8 * maxfact))), sep = "", "\n")#ラインを引く；要因の数に合わせて長さを調整
	cat(sprintf("%3s", names(bstatist[1:maxfact])), sprintf("%2s", names(bstatist[maxfact+1])), sprintf("%6s", names(bstatist[maxfact+2])), 
		sprintf("%14s", names(bstatist[maxfact+3])), sep = "\t", "\n")#データフレームの列名をタブ区切りで表示する
	cat(sprintf(rep("-", (36 + 8 * maxfact))), sep = "", "\n")

	#データフレームのデータ部分を一行ずつタブ区切りで表示する
	for (i in 1:nrow(bstatist)){
		if((i %% margin == 0) && !(i == nrow(bstatist))){
			cat(sprintf("%3s", unlist(bstatist[i, 1:maxfact])), sprintf("%3s", paste(bstatist[i, maxfact+1])), sprintf("%9.4f", bstatist[i, maxfact+2:3]), sep = "\t", "\n")
			cat("\n")
		}else cat(sprintf("%3s", unlist(bstatist[i, 1:maxfact])), sprintf("%3s", paste(bstatist[i, maxfact+1])), sprintf("%9.4f", bstatist[i, maxfact+2:3]), sep = "\t", "\n")
	}

	cat(sprintf(rep("-", (36 + 8 * maxfact))), sep = "", "\n", "\n", "\n")

}


#分散分析表を出力する関数
anovatab.out <- function(mainresults){
	#情報を展開
	epsi.info1 <- mainresults$epsi.info1
	epsitab <- mainresults$epsitab
	ano.info1 <- mainresults$ano.info1
	anovatab <- mainresults$anovatab
	mse.row <- mainresults$mse.row

	#球面性検定の結果を出力
	if(is.na(epsi.info1) == FALSE){#epsi.info1がNAでないときのみ出力
		cat(sprintf("<< SPHERICITY INDICES >>"), sep = "", "\n", "\n")#タイトル
		cat(sprintf(epsi.info1), "\n", "\n")#プロンプト

		cat(sprintf(rep("-", 6 + 9 * ncol(epsitab))), sep = "", "\n")#ラインを引く
		cat(sprintf("%10s", names(epsitab[1])), sprintf("%6s", names(epsitab[2])), sprintf("%8s", names(epsitab[3])), sprintf("%3s", names(epsitab[4])), 
			sprintf("%3s", names(epsitab[5])), sprintf(""), sprintf("%3s", names(epsitab[7:ncol(epsitab)])), sep = "\t", "\n")
		cat(sprintf(rep("-", 6 + 9 * ncol(epsitab))), sep = "", "\n")

		#データフレームのデータ部分を一行ずつタブ区切りで表示する
		for (i in 1:nrow(epsitab)){
			cat(sprintf("%10s", paste(epsitab[i, 1])), sprintf("%7.4f", epsitab[i, 2]), sprintf("%8.4f", epsitab[i, 3]), 
				sprintf("%3.0f", epsitab[i, 4]), replace(sprintf("%6.4f", epsitab[i, 5]), is.na(epsitab[i, 5]), ""), 
				sprintf("%4s", paste(epsitab[i, 6])), sprintf("%6.4f", epsitab[i, 7:ncol(epsitab)]), sep = "\t", "\n")
		}

		cat(sprintf(rep("-", 6 + 9 * ncol(epsitab))), sep = "", "\n")#ラインを引く
		if(ncol(epsitab) == 9){
			cat(sprintf(paste(rep("", 3), sep = "")), sprintf("%62s", "LB = lower.bound, GG = Greenhouse-Geisser, HF = Huynh-Feldt"), sep = "\t", "\n", "\n", "\n")
		}else{
			cat(sprintf(paste(rep("", 3), sep = "")), sprintf("%78s", "b_hat = b^, c_hat = c^, h_d = h~', h_dd = h~'', h = h~"), sep = "\t", "\n", "\n", "\n")
		}
	}

	#分散分析表を出力
	cat(sprintf("<< ANOVA TABLE >>"), sep = "", "\n", "\n")#タイトル
	if(sum(is.na(ano.info1)) != 1){#ano.info1がNAでないときのみ表示
		cat(sprintf(ano.info1), sep = "\n")#プロンプトを表示
		cat("\n")#１行空ける
	}

	if(ncol(anovatab) == 7){#分散分析表のみの場合
		cat(sprintf(rep("-", 84)), sep = "", "\n")#ラインを引く
		cat(sprintf("%10s", "Source"), sprintf("%16s", paste("SS")), sprintf("%12s", "df"), sprintf("%14s", "MS"), 
			sprintf("%12s", "F-ratio"), sprintf("%16s", paste("p-value", "", sep = "\t")), sep = "", "\n")#列名をタブ区切りで表示する
		cat(sprintf(rep("-", 84)), sep = "", "\n")

		#データフレームのデータ部分を一行ずつタブ区切りで表示する；誤差項の行（mse.row）の後にラインを入れる
		for (i in 1:nrow(anovatab)){
			if(i == nrow(anovatab)){#最終行に到達したときには以下を実行
				cat(sprintf("%10s", paste(anovatab[i, 1])), sprintf("%14.4f", anovatab[i, 2]), 
					sprintf("%6s", round(anovatab[i, 3], 2)), sprintf("%-8s", "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = "\t", "\n")
			}else if(sum(mse.row == i) == 1){#mse.rowの中に一致する値があるときには以下を実行
				cat(sprintf("%10s", paste(anovatab[i, 1])), sprintf("%14.4f", anovatab[i, 2]), 
					sprintf("%6s", round(anovatab[i, 3], 2)), sprintf("%14.4f", anovatab[i, 4]), 
					replace(sprintf("%8.4f", anovatab[i, 5]), is.na(anovatab[i, 5]), ""), 
					replace(sprintf("%.4f", anovatab[i, 6]), is.na(anovatab[i, 6]), ""), sprintf(paste(anovatab[i, 7])), sep = "\t", "\n")
				cat(sprintf(rep("-", 84)), sep = "", "\n")
			}else{#その他のときは以下を実行
				cat(sprintf("%10s", paste(anovatab[i, 1])), sprintf("%14.4f", anovatab[i, 2]), 
					sprintf("%6s", round(anovatab[i, 3], 2)), sprintf("%14.4f", anovatab[i, 4]), 
					replace(sprintf("%8.4f", anovatab[i, 5]), is.na(anovatab[i, 5]), ""), replace(sprintf("%.4f", anovatab[i, 6]), 
					is.na(anovatab[i, 6]), ""), sprintf(paste(anovatab[i, 7])), sep = "\t", "\n")
			}
		}
	}else{#効果量の指標を追加した場合
		esn <- ncol(anovatab) - 7
		cat(sprintf(rep("-", 84 + 16 * esn)), sep = "", "\n")#ラインを引く
		cat(sprintf("%10s", "Source"), sprintf("%16s", paste("SS")), sprintf("%12s", "df"), sprintf("%14s", "MS"), 
			sprintf("%12s", "F-ratio"), sprintf("%16s", paste("p-value", "", sep = "\t")), sprintf("%16s", names(anovatab[8:ncol(anovatab)])), sep = "", "\n")#列名をタブ区切りで表示する
		cat(sprintf(rep("-", 84 + 16 * esn)), sep = "", "\n")

		#データフレームのデータ部分を一行ずつタブ区切りで表示する；誤差項の行（mse.row）の後にラインを入れる
		for (i in 1:nrow(anovatab)){
			if(i == nrow(anovatab)){#最終行に到達したときには以下を実行
				cat(sprintf("%10s", paste(anovatab[i, 1])), sprintf("%14.4f", paste(anovatab[i, 2])), 
					sprintf("%6s", round(anovatab[i, 3], 2)), sprintf(paste("%", 40 + 16 * esn, "s", sep = ""), "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = "\t", "\n")
			}else if(sum(mse.row == i) == 1){#mse.rowの中に一致する値があるときには以下を実行
				cat(sprintf("%10s", paste(anovatab[i, 1])), sprintf("%14.4f", anovatab[i, 2]), 
					sprintf("%6s", round(anovatab[i, 3], 2)), sprintf("%14.4f", anovatab[i, 4]), 
					replace(sprintf("%8.4f", anovatab[i, 5]), is.na(anovatab[i, 5]), ""), 
					replace(sprintf("%.4f", anovatab[i, 6]), is.na(anovatab[i, 6]), ""), sprintf(paste(anovatab[i, 7])), sep = "\t", "\n")
				cat(sprintf(rep("-", 84 + 16 * esn)), sep = "", "\n")
			}else{#その他のときは以下を実行
				cat(sprintf("%10s", paste(anovatab[i, 1])), sprintf("%14.4f", anovatab[i, 2]), 
					sprintf("%6s", round(anovatab[i, 3], 2)), sprintf("%14.4f", anovatab[i, 4]), 
					replace(sprintf("%8.4f", anovatab[i, 5]), is.na(anovatab[i, 5]), ""), replace(sprintf("%.4f", anovatab[i, 6]), 
					is.na(anovatab[i, 6]), ""), sprintf(paste(anovatab[i, 7])), replace(sprintf("%8.4f", anovatab[i, 8:ncol(anovatab)]), is.na(anovatab[i, 8:ncol(anovatab)]), ""), sep = "\t", "\n")
			}
		}
	}

	cat("\n", "\n")#２行空ける

}


#修正Bonferroniの方法による多重比較の結果を出力する関数
mod.Bon.out <- function(bon.list, omit = FALSE){
	factlabel <- bon.list[[1]]
	bon.info1 <- bon.list[[2]]
	bon.info2 <- bon.list[[3]]
	bon.info3 <- bon.list[[4]]
	bonstat <- bon.list[[5]]
	bontab <- bon.list[[6]]

	cat(sprintf(paste("< MULTIPLE COMPARISON for FACTOR ", factlabel, " >" ,sep = "")), sep = "", "\n", "\n")#タイトル
	cat(sprintf(bon.info1), sep = "\n")#プロンプト
	cat(sprintf(bon.info2), "\n")
	cat(sprintf(bon.info3), "\n", "\n")

	#omitがFなら記述統計量を出力する
	if(omit == FALSE){
		#記述統計量を出力する
		cat(sprintf(rep("-", 44)), sep = "", "\n")#ラインを引く；要因の数に合わせて長さを調整
		cat(sprintf("%4s", names(bonstat)[1]), sprintf("%2s", names(bonstat)[2]), sprintf("%6s", names(bonstat)[3]), 
			sprintf("%14s", names(bonstat)[4]), sep = "\t", "\n")#データフレームの列名をタブ区切りで表示する
		cat(sprintf(rep("-", 44)), sep = "", "\n")

		#データフレームのデータ部分を一行ずつタブ区切りで表示する
		for (i in 1:nrow(bonstat)){
			cat(sprintf("%4s", paste(bonstat[i, 1])), sprintf("%2s" ,paste(bonstat[i, 2])), sprintf("%9.4f", bonstat[i, 3]), 
				sprintf("%9.4f", bonstat[i, 4]), sep = "\t", "\n")
		}

		cat(sprintf(rep("-", 44)), sep = "", "\n", "\n")	

	}	

	#多重比較の結果を表形式で出力する
	cat(sprintf(rep("-", 75)), sep = "", "\n")
	cat(sprintf("%6s", "Pair"), sprintf("%12s", paste("Interval")), sprintf("%14s", "t-value"), sprintf("%12s", "df"), 
		sprintf("%7s", "p"), sprintf(paste("%", nchar(names(bontab[6]))+5, "s"), names(bontab[6])), sep = "", "\n")#列名をタブ区切りで表示する
	cat(sprintf(rep("-", 75)), sep = "", "\n")

	#データフレームのデータ部分を一行ずつタブ区切りで表示する
	for (i in 1:nrow(bontab)){
		cat(sprintf("%6s",paste(bontab[i, 1])), sprintf("%10.4f", bontab[i, 2]), sprintf("%8.4f", bontab[i, 3]), 
			sprintf(paste(" ", bontab[i, 4])), sprintf("%6.4f", bontab[i, 5]), sprintf("%6.4f", bontab[i, 6]), 
			sprintf(paste(bontab[i, 7])), sep = "\t", "\n")
	}

	cat(sprintf(rep("-", 75)), sep = "", "\n", "\n")

}


#単純主効果の検定の結果を出力する関数
simeffects.out <- function(partresult, omit = FALSE){
	part.info1 <- partresult$sig.term
	partstat <- partresult$sim.stat
	partepsi <- partresult$simepsi
	parttab <- partresult$simtab
	partmulttab <- partresult$sim.multresults

	cat(sprintf(paste("< SIMPLE EFFECTS for <", part.info1, "> INTERACTION >"), sep = ""), sep = "", "\n", "\n")#タイトル

	#omitがFALSEなら記述統計量を出力する
	if(omit == FALSE){
		#記述統計量を出力する
		cat(sprintf(rep("-", 52)), sep = "", "\n")#ラインを引く；要因の数に合わせて長さを調整
		cat(sprintf("%4s", names(partstat)[1]), sprintf("%2s", names(partstat)[2]), sprintf("%s", names(partstat)[3]), 
			sprintf("%6s", names(partstat)[4]), sprintf("%14s", names(partstat)[5]), sep = "\t", "\n")#データフレームの列名をタブ区切りで表示する
		cat(sprintf(rep("-", 52)), sep = "", "\n")

		#データフレームのデータ部分を一行ずつタブ区切りで表示する
		for (i in 1:nrow(partstat)){
			cat(sprintf("%4s", paste(partstat[i, 1])), sprintf(paste(partstat[i, 2])), sprintf(paste(partstat[i, 3])), sprintf("%9.4f", partstat[i, 4:5]), sep = "\t", "\n")
		}

		cat(sprintf(rep("-", 52)), sep = "", "\n", "\n")	

	}

	#球面性検定の結果を出力
	if(is.na(charmatch("Error", parttab$source.col)) && !is.na(partepsi)){
		#被験者内計画，混合要因計画なら球面性検定の結果を出力する
		cat(sprintf(rep("-", 6 + 9 * ncol(partepsi))), sep = "", "\n")#ラインを引く
		cat(sprintf("%10s", names(partepsi[1])), sprintf("%6s", names(partepsi[2])), sprintf("%8s", names(partepsi[3])), sprintf("%3s", names(partepsi[4])), 
			sprintf("%3s", names(partepsi[5])), sprintf(""), sprintf("%5s", names(partepsi[7:ncol(partepsi)])), sep = "\t", "\n")
		cat(sprintf(rep("-", 6 + 9 * ncol(partepsi))), sep = "", "\n")

		#データフレームのデータ部分を一行ずつタブ区切りで表示する
		for (i in 1:nrow(partepsi)){
			cat(sprintf("%10s", paste(partepsi[i, 1])), sprintf("%7.4f", partepsi[i, 2]), 
				sprintf("%8.4f", partepsi[i, 3]), sprintf("%3.0f", partepsi[i, 4]), 
				replace(sprintf("%6.4f", partepsi[i, 5]), is.na(partepsi[i, 5]), ""), 
				sprintf("%4s", paste(partepsi[i, 6])), sprintf("%6.4f", partepsi[i, 7:ncol(partepsi)]), sep = "\t", "\n")
		}

		cat(sprintf(rep("-", 6 + 9 * ncol(partepsi))), sep = "", "\n")#ラインを引く
		if(ncol(partepsi) == 9){
			cat(sprintf(paste(rep("", 3), sep = "")), sprintf("%62s", "LB = lower.bound, GG = Greenhouse-Geisser, HF = Huynh-Feldt"), sep = "\t", "\n", "\n")
		}else{
			cat(sprintf(paste(rep("", 3), sep = "")), sprintf("%78s", "b_hat = b^, c_hat = c^, h_d = h~', h_dd = h~'', h = h~"), sep = "\t", "\n", "\n")
		}

	}

	if(ncol(parttab) == 7){#分散分析表のみの場合
		#単純主効果の検定の分散分析表を出力する
		cat(sprintf(rep("-", 84)), sep = "", "\n")#ラインを引く
		cat(sprintf("%10s", "Source"), sprintf("%16s", paste("SS")), sprintf("%12s", "df"), sprintf("%14s", "MS"), 
			sprintf("%12s", "F-ratio"), sprintf("%16s", paste("p-value", "", sep = "\t")), sep = "", "\n")#列名をタブ区切りで表示する
		cat(sprintf(rep("-", 84)), sep = "", "\n")

		#データフレームのデータ部分を一行ずつタブ区切りで表示する
		if(is.na(charmatch("Error", parttab$source.col))){
			#被験者内計画，混合要因計画なら２行ごとにラインを入れる
			for (i in 1:nrow(parttab)){
				cat(sprintf("%10s", paste(parttab[i, 1])), sprintf("%14.4f", parttab[i, 2]), 
					sprintf("%6s", round(parttab[i, 3], 2)), sprintf("%14.4f", parttab[i, 4]), 
					replace(sprintf("%8.4f", parttab[i, 5]), is.na(parttab[i, 5]), ""), 
					replace(sprintf("%.4f", parttab[i, 6]), is.na(parttab[i, 6]), ""), sprintf(paste(parttab[i, 7])), sep = "\t", "\n")
				if((i %% 2) == 0){#行番号が２で割り切れるときはラインを引く
					cat(sprintf(rep("-", 84)), sep ="", "\n")
				}
			}
			#表の末尾に有意性判定に使った記号を記す
			cat(sprintf(paste(rep("", 5), sep = "")), sprintf("%-8s", "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = "\t", "\n", "\n")
		}else{
			#被験者間計画なら誤差項の上の行だけにラインを入れる
			for (i in 1:nrow(parttab)){
				if(nrow(parttab) == i){#行番号が最後の行に一致するときは以下を実行
					cat(sprintf("%10s", paste(parttab[i, 1])), sprintf("%14.4f", paste(parttab[i, 2])), sprintf("%6s", paste(parttab[i, 3])), 
						sprintf("%14.4f", paste(parttab[i, 4])), sep = "\t", "\n")
					cat(sprintf(rep("-", 84)), sep = "", "\n")
					cat(sprintf(paste(rep("", 5), sep = "")), sprintf("%-8s", "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = "\t", "\n", "\n")
				}else{#その他のときは以下を実行
					cat(sprintf("%10s", paste(parttab[i, 1])), sprintf("%14.4f", parttab[i, 2]), sprintf("%6s", paste(parttab[i, 3])), 
						sprintf("%14.4f", parttab[i, 4]), replace(sprintf("%8.4f", parttab[i, 5]), is.na(parttab[i, 5]), ""), 
						replace(sprintf("%.4f", parttab[i, 6]), is.na(parttab[i, 6]), ""), sprintf(paste(parttab[i, 7])), sep = "\t", "\n")
				}
			}
		}
	}else{#効果量の指標を追加した場合
		esn <- ncol(parttab) - 7
		#単純主効果の検定の分散分析表を出力する
		cat(sprintf(rep("-", 84 + 16 * esn)), sep = "", "\n")#ラインを引く
		cat(sprintf("%10s", "Source"), sprintf("%16s", paste("SS")), sprintf("%12s", "df"), sprintf("%14s", "MS"), 
			sprintf("%12s", "F-ratio"), sprintf("%16s", paste("p-value", "", sep = "\t")), sprintf("%16s", names(parttab[8:ncol(parttab)])), sep = "", "\n")#列名をタブ区切りで表示する
		cat(sprintf(rep("-", 84 + 16 * esn)), sep = "", "\n")

		#データフレームのデータ部分を一行ずつタブ区切りで表示する
		if(is.na(charmatch("Error", parttab$source.col))){
			#被験者内計画，混合要因計画なら２行ごとにラインを入れる
			for (i in 1:nrow(parttab)){
				cat(sprintf("%10s", paste(parttab[i, 1])), sprintf("%14.4f", parttab[i, 2]), 
					sprintf("%6s", round(parttab[i, 3], 2)), sprintf("%14.4f", parttab[i, 4]), 
					replace(sprintf("%8.4f", parttab[i, 5]), is.na(parttab[i, 5]), ""), 
					replace(sprintf("%.4f", parttab[i, 6]), is.na(parttab[i, 6]), ""), sprintf(paste(parttab[i, 7])), 
					replace(sprintf("%8.4f", parttab[i, 8:ncol(parttab)]), is.na(parttab[i, 8:ncol(parttab)]), ""), sep = "\t", "\n")
				if((i %% 2) == 0){#行番号が２で割り切れるときはラインを引く
					cat(sprintf(rep("-", 84 + 16 * esn)), sep ="", "\n")
				}
			}
			#表の末尾に有意性判定に使った記号を記す
			cat(sprintf(paste(rep("", 5), sep = "")), sprintf(paste("%", 40 + 16 * esn, "s", sep = ""), "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = "\t", "\n", "\n")
		}else{
			#被験者間計画なら誤差項の上の行だけにラインを入れる
			for (i in 1:nrow(parttab)){
				if(nrow(parttab) == i){#行番号が最後の行に一致するときは以下を実行
					cat(sprintf("%10s", paste(parttab[i, 1])), sprintf("%14.4f", paste(parttab[i, 2])), sprintf("%6s", paste(parttab[i, 3])), 
						sprintf("%14.4f", paste(parttab[i, 4])), sep = "\t", "\n")
					cat(sprintf(rep("-", 84 + 16 * esn)), sep = "", "\n")
					cat(sprintf(paste(rep("", 5), sep = "")), sprintf(paste("%", 40 + 16 * esn, "s", sep = ""), "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = "\t", "\n", "\n")
				}else{#その他のときは以下を実行
					cat(sprintf("%10s", paste(parttab[i, 1])), sprintf("%14.4f", parttab[i, 2]), sprintf("%6s", paste(parttab[i, 3])), 
						sprintf("%14.4f", parttab[i, 4]), replace(sprintf("%8.4f", parttab[i, 5]), is.na(parttab[i, 5]), ""), 
						replace(sprintf("%.4f", parttab[i, 6]), is.na(parttab[i, 6]), ""), sprintf(paste(parttab[i, 7])), 
						replace(sprintf("%8.4f", parttab[i, 8:ncol(parttab)]), is.na(parttab[i, 8:ncol(parttab)]), ""), sep = "\t", "\n")
				}
			}
		}

	}

	#多重比較の結果を出力する
	for (i in 1:length(partmulttab)){
		if(is.na(partmulttab[[i]][1]) == FALSE){#リストに中身があったら出力する
			cat("\n")#１行空ける
			mod.Bon.out(partmulttab[[i]], omit = TRUE)
		}
	}

	cat("\n")#１行空ける

}


#下位検定の結果を出力する関数
post.out <- function(postresults, design){
	#リストに含まれる要素の数を特定
	postnum <- length(postresults)

	#すべてのリストの中身がNAの場合には，タイトルも表示しない
	if(sum(is.na(postresults)) != postnum){
		cat(sprintf("<< POST ANALYSES >>"), sep = "", "\n", "\n")#タイトル
	}

	#リストの中身をひとつずつeach.out関数に送る
	for (i in 1:postnum){
		each.out(postresults[[i]], design)
	}

	#出力が終わったことを知らせるプロンプト
	cat(sprintf("output is over "), sprintf(rep("-", 20)), sprintf("///"), sep = "", "\n")
	cat("\n")#１行空ける
}


#効果のタイプによって異なる出力方式を割り当てる関数
each.out <- function(partresult, design){
	if(is.na(partresult[[1]]) == TRUE){
	#何もない場合（２水準の主効果が有意で多重比較の必要なしの場合）

	}else if(nchar(partresult[[1]]) == 1){
	#多重比較の結果がある場合
		mod.Bon.out(partresult)
		cat("\n")#１行空ける

	}else if(nchar(partresult[[1]]) == 3){
	#一次の交互作用が見られた場合
		if(nchar(design) == 3) simeffects.out(partresult, omit = TRUE)#もともと２要因の計画のときには，単純主効果の検定の出力時に記述統計量を出力しない（主分析の記述統計量と重複するので）
		else simeffects.out(partresult)

	}else if(nchar(partresult[[1]]) >= 5){
	#高次の交互作用が見られた場合
		cat(sprintf(paste("< HIGHER-ORDER < ", partresult[[1]], " > INTERACTION >", sep = "")),sep = "", "\n")
		cat(sprintf("*** Split the dataset for further analysis. ***"), sep ="", "\n", "\n")#データ分割を促すプロンプト
	}
}

