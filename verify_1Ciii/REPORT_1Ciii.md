# 1-C-iii g-forcing 検証レポート（nm=32 段階検証）

**ジョブ**: 2353404（Miyabi `regular-g`, 16 GPU, `P_fft=4 / P_m=4 / P_s=1`）
**解像度**: 64² × 32, Hermite `nm=32`
**日付**: 2026-07-10

## 目的

prod128 の**フリーディケイ問題**（Alfvén のみ強制で圧縮性 g カスケードに自由エネルギー源が無く、
移流 `{Φ,g}` と平行ストリーミングが共に `W=½Σ|g_m|²` を保存 → W が単調減衰するのみ）を、
Meyrand 2019 に倣い **g を m=1・低波数で一定の圧縮性自由エネルギー注入率で駆動**する修正で解消する。
本 run はその修正を**本番 128²×64 の前に安価に検証**する段階ゲート。

## 使用ファイル

| 種別 | ファイル | 役割 |
|------|----------|------|
| 入力 (driven) | `echo64.in` | nonlinear=**T** 駆動乱流。フリーディケイ修正の本番。 |
| 入力 (control) | `lin64.in` | nonlinear=**F** 線形コントロール（前方 Hermite カスケード）。`echo64.in` と max_wall_time 以外同一。 |
| ジョブスクリプト | `job.pbs-1Ciii` | 同一バイナリで lin64 → echo64 を連続実行（lin64 が NetCDF を出さなければ echo64 はスキップ）。 |
| 解析 (remote) | `analyze.py` | ジョブ内で ncdump 経由の数値サマリを出力。 |
| 作図 (local) | `plot_1Ciii.py` | 本レポートの図を NetCDF から local 生成。 |

### 主要パラメータ（`echo64.in`）

- `&force`: `driven=T, elsasser=.true., fix_power=.true., ene_inj=1.0, xhl_inj=0.8, **ene_inj_g=1.0**, kmin=(0,0,1), kmax=(1,1,2), nfields=3, stir_seed=7`
- `&forced_fields`: `field_names='zppe','zmpe','g1'`（**g1 を第3の強制場として Alfvén と同じ低 k 撹拌シェルで駆動**）, `frequencies=(0.9,-0.6)` ×3
- `&physical_parameters`: `v_th=1, mu_hyper_perp=1d1 (nexp_perp=4), nu_hyper_m=1d0 (nexp_m=6), beta_i=tau=Zcharge=1, alpha_root=1`（散逸係数は prod128 と同一・k_max/m=nm 正規化で解像度非依存）
- `&time_parameters`: `dt=2e-4, cfl=0.5, reset_method='decrement'`
- `&operation_parameters`: `nonlinear=T, write_hermite_flux=T`

## 結果

### 図1: 駆動定常 W_free(t)

![W_free(t)](figs/fig_W_free.png)

両 run とも `W_free` が初期 5.44 から**プラトーへ上昇し揺動**する（lin64: 平均 9.89, echo64: 平均 9.12）。
prod128 で見られた **max=初期値の単調減衰は消滅**しており、g 強制が圧縮性自由エネルギーを正しく注入している。

### 図2: Hermite 自由エネルギースペクトル W_m(m)

![W_m](figs/fig_W_m.png)

両者とも強制モード m=1 でピーク後に減衰。**非線形 echo64（青）は線形 lin64（橙）より急峻**で、
中域 m≈4–15 で `m⁻¹`（fluidization）側に寄り、lin64 は `m⁻¹/²`（純位相混合）側に沿う。
高 m（>15）は偶奇パリティ振動が支配し nm=32 では未分解。

### 図3: Hermite フラックス Γ_m と echo metric

![Gamma_m](figs/fig_Gamma_m.png)

- **左（k⊥=7 の累積 Γ_m）**: lin64（橙）はほぼゼロ — 非線形が無いと**perp カスケードが起きず g エネルギーは強制した低 k⊥ に留まる**ため k⊥=7 は空。echo64（青）は低 m で前方（正）→ 中域で単調減衰（anti-mixing）→ 高 m で負（散逸域）。
- **右（inertial band の ECHO_RATIO）**: lin64 ≈ 0.44 でほぼ平坦、echo64 ≈ 0.58。nm=32 では両者とも O(0.5) で、clean な echo（<<1）には未到達。

### 数値サマリ（`analyze.py`, steady window t≥4.0）

| 指標 | lin64 (nonlinear=F) | echo64 (nonlinear=T) | 判定 |
|------|--------------------|--------------------|------|
| NREC / T_REACHED | 634 / 158.1 | 98 / 24.2 | — |
| W_FREE_FIRST → MEAN | 5.44 → **9.89** | 5.44 → **9.12** | ✅ プラトー |
| W_FREE_STEADY_RELRMS | **0.0426** | **0.0750** | ✅ 揺動（単調減衰でない） |
| TELE_KINT (Σ src_m) | 6.71e-16 | 6.26e-16 | ✅ ≈0（telescoping 健全） |
| ECHO_RATIO_MEAN | 0.4396 | 0.5767 | ⚠️ nm=32 では echo 未分解 |

## 合否判定

| ゲート | 内容 | 結果 |
|--------|------|------|
| **(1) 駆動定常＝フリーディケイ修正** | W_free がプラトー＋揺動、単調減衰でない | ✅ **合格**（本 run の主目的） |
| **telescoping sanity** | Σ_m src_m ≈ 0（Γ_m 診断の自己整合） | ✅ **合格**（~6e-16, 機械精度） |
| 前方カスケード contrast | lin64 で forward flux（echo 無し） | ✅ 確認 |
| echo 方向性のヒント | echo64 の W_m が急峻・Γ_m が中域で減衰 | 🟡 定性的に echo 方向（決定的でない） |
| **(2) stochastic echo（Γ_m≈0, m⁻¹）** | inertial band で ECHO_RATIO<<1 | ⏳ **未分解 → 本番 nm=128 待ち** |

## 結論

**フリーディケイ問題は解消**した（g 強制の配線・注入率正規化・3 段 RK 全ステージ注入が正しく機能し、
両 run が駆動定常に到達）。telescoping sanity も機械精度で合格。段階計画通り、
**漸近的な stochastic echo（Γ_m≈0 の inertial band・m⁻¹ スペクトル）は Hermite 慣性領域が短い nm=32 では未分解**で、
本番 **128²×64, nm=128** の高解像度 run が次段。echo 方向のヒント（echo64 の W_m の急峻化・Γ_m の中域減衰）は励みになる。

**次段**: prod128 の入力を driven-g 仕様（`ene_inj_g`, `field_names` に `g1` 追加）に更新し、128²×64/nm=128 で本番 echo 図を取得。
