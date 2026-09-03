# physmerge — 技術規格 (Technical Specification)

- **Package**: physmerge
- **Version**: 0.3.0
- **License**: MIT
- **Language**: R (≥ base)
- **Imports**: `utils` (zip)
- **Suggests**: `data.table` (檔案讀取)、`shiny`/`DT`/`shinyjs` (未來 UI)、`testthat (≥ 3.0.0)`
- **Author**: Pin
- **Repo**: `Droideight/physmerge`

---

## 1. 套件目標 (Goal)

將 GWAS summary statistics 中**位置上相鄰**的顯著訊號 (significant SNPs)，以一個**前向滑動視窗 (forward sliding-window)** 演算法塌縮成**互不重疊** (zero-overlap) 的 locus blocks。

設計核心：
- **不需要 LD reference panel** — 只用 base-pair 距離判斷鄰近性
- **單次線性掃描** (O(n) 在 sort 之後)，記憶體連續配置
- **per-chromosome 自動分流**，防止染色體邊界錯誤合併
- **可同時處理 p-value (越小越顯著) 與 test statistic (越大越顯著)**，透過 `reward ∈ {"min","max"}` 抽象化

---

## 2. 公開 API (Exported Functions)

| 函式 | 檔案 | 簽章 |
|---|---|---|
| `read_sumstat()` | [R/read_sumstat.R](R/read_sumstat.R) | 讀檔 + 欄位正規化 → `list(data, reward)` |
| `physical_merge()` | [R/physical_merge.R](R/physical_merge.R) | 合併核心 → blocks data frame |
| `annotate_blocks()` | [R/export.R](R/export.R) | 把原始欄位接回 blocks |
| `export_snp_list()` | [R/export.R](R/export.R) | 輸出代表性 SNP 清單 (`.txt` 或 per-chrom `.zip`) |

另有一支獨立的 C99 執行檔 [cli/physmerge.c](cli/physmerge.c)，把上述四個函式合成單一指令，不需要 R。串流單次掃描，記憶體固定 2.4 MB。每個 R 參數都有對應 flag，對照表與已知差異見 [cli/README.md](cli/README.md)、效能與驗證見 [cli/PERFORMANCE.md](cli/PERFORMANCE.md)。

內部 helpers (非 exported)：
- `.is_significant(val, sig_th, reward)` — 顯著性比較
- `.is_more_significant(val, best, reward)` — 「更顯著」比較
- `.physical_merge_single(...)` — 單染色體核心掃描
- `.collapse_blocks(blk, w, reward)` — 鄰近 blocks 合併
- `.stash_rps_row(blk)` — 把代表 SNP 的輸入列號搬到 attribute
- `.finish_annotation(...)` — `annotate_blocks()` 兩條路徑共用的欄位排序
- `` `%||%` `` — null-coalescing operator

---

## 3. 資料契約 (Data Contract)

### 3.1 輸入：標準化欄位
所有合併皆基於兩個必要欄位：
- `position` (numeric, bp 座標)
- `value` (numeric, p-value 或 test statistic)

可選欄位：
- `CHROM` (character/integer/numeric) — 觸發 per-chromosome 分流
- `ID` / `SNP` — 代表性 SNP 識別碼，供 `annotate_blocks()` 使用

### 3.2 輸出：Blocks data frame

| 欄位 | 型別 | 說明 |
|---|---|---|
| `serial` | integer | 連續 block 編號 (跨染色體一次重編) |
| `CHROM` | (取自原始) | 僅在偵測到染色體欄位時出現 |
| `start` | numeric | block 起始 bp，= `max(0, rps_BP_first − window)` |
| `end` | numeric | block 結束 bp，= 最後 in-block SNP 位置 + 剩餘 `steps` |
| `rps_BP` | numeric | Representative SNP 位置 (最顯著 SNP 的 bp) |
| `rps_value` | numeric | Representative SNP 的 `value` |

此外攜帶一個 **attribute** `attr(blocks, "rps_row")`：每個 block 的代表 SNP 在輸入 `data` 中的**列號** (integer vector，長度等於 `nrow(blocks)`)。它不是欄位，不影響輸出長相；`annotate_blocks()` 用它精準取回代表列，解決同 bp 多 SNP 的標註問題 (見 6.2)。若使用者對 blocks 做了列的增刪，長度或位置檢核會失敗，`annotate_blocks()` 自動退回舊的座標比對。

不變式 (invariants)：
- `start[i] ≤ rps_BP[i] ≤ end[i]`
- 對任意連續兩個 blocks：`end[i] ≤ start[i+1]` (Trim pass 保證)
- 輸出依 `position` 升冪排序 (per chromosome)

---

## 4. 核心演算法 — `physical_merge()`

### 4.1 三階段管線 (Three-pass pipeline)

```
       ┌──────────────────────────────────────────┐
input  │ Pass 1: Forward scan (.physical_merge_single)
data ──▶│   sort by position                     │
       │   open/extend/close raw blocks         │
       └──────────────┬───────────────────────────┘
                      ▼
       ┌──────────────────────────────────────────┐
       │ Pass 2: Collapse (.collapse_blocks)     │
       │   merge if Δrps_BP < window             │
       └──────────────┬───────────────────────────┘
                      ▼
       ┌──────────────────────────────────────────┐
       │ Pass 3: Trim (inline in single)         │
       │   if end[i] > start[i+1]: end[i]=start[i+1]
       └──────────────┬───────────────────────────┘
                      ▼
                 final blocks
```

### 4.2 Pass 1 — Forward scan (狀態機)

**狀態變數**：
- `in_block` (logical) — 目前是否在 block 內
- `steps` (numeric) — 視窗剩餘餘額 (bp)
- `sig_this_block` — 目前 block 內最顯著的 value
- `last_pos` — 前一筆 SNP 位置

**單筆迭代邏輯** (`for i in 1..n`)：

```
case (!in_block):
    if value 顯著:  open_block(pos, val)

case (in_block):
    remaining = steps − (pos − last_pos)

    case remaining <= 0:                  # 視窗耗盡
        close_block(last_pos)
        if value 顯著: open_block(pos, val)

    case remaining > 0:                    # 仍在視窗內
        steps = remaining
        if reset_on == "any" 且 value 顯著:
            steps = window                 # 任何顯著 SNP 都重置
            if value 更顯著:
                update representative
        else if value 更顯著:               # reset_on == "best"
            steps = window
            update representative
```

收尾：若迴圈結束時仍 `in_block`，呼叫 `close_block(last_pos)`。

**`open_block(pos, val)`**：
- `block_count += 1`
- `start = max(0, pos − window)` ← 向左外推一個 window
- `rps_BP = pos`, `rps_value = val`
- 重置 `steps = window`、`sig_this_block = val`

**`close_block(last_inblock_pos)`**：
- `end = last_inblock_pos + steps` ← 向右延伸剩餘餘額
- 重置狀態

### 4.3 `reset_on` 兩種語意 (semantics)

| 模式 | 觸發視窗重置的條件 | 等價於 |
|---|---|---|
| `"best"` (預設) | 只在遇到**更顯著**的 SNP 時重置 | 原 physmerge 行為，代表 SNP 永遠是 local optimum |
| `"any"` | 任何**顯著** (≥ `sig_th`) SNP 都重置 | locusDefiner：取所有顯著 SNP 的 ±window 區間之聯集 |

注意：`"any"` 模式下，代表 SNP 仍只在遇到「**更**顯著」時才更新；window 重置與代表更新是**分離**的兩件事。

### 4.4 Pass 2 — Collapse pass (`.collapse_blocks`)

線性掃描已產生的 raw blocks：

```
out = [blk[1]]
for i in 2..nrow(blk):
    if (blk[i].rps_BP − out[last].rps_BP) < window:
        out[last].end = max(out[last].end, blk[i].end)
        if blk[i].rps_value 更顯著:
            out[last].rps_BP    = blk[i].rps_BP
            out[last].rps_value = blk[i].rps_value
    else:
        out.append(blk[i])
重編 serial
```

合併條件**依代表 SNP 距離** (不是 block 端點距離)，因此可吃下「主峰被次峰拉開」的情況。

### 4.5 Pass 3 — Trim pass

```
for i in 1..(nrow(blk) − 1):
    if blk$end[i] > blk$start[i+1]:
        blk$end[i] = blk$start[i+1]
```

確保任何重疊都被消除 (zero-overlap invariant)。

### 4.6 Per-chromosome 分流

於 `physical_merge()` (top-level wrapper)：

1. **解析染色體欄名**：使用者 `chrom_col` 參數 → fallback 到 `"CHROM"` → 否則 `NULL`
2. 若 `length(unique(data[[chrom_col]])) > 1`：
   - `lapply(chroms, ...)` 分群執行 `.physical_merge_single()`
   - `rbind` 後重編 `serial`
   - 把 `CHROM` 欄位插到 `serial` 之後
3. 若**無**染色體欄位但 `range(position) > 250 Mb`：
   - 發出 `warning()`，提醒可能跨染色體誤合併

### 4.7 輸入驗證 (Input validation)

於 `physical_merge()` 進入點：
- `is.data.frame(data)`
- `c("position","value") %in% names(data)`
- `is.numeric` for both
- `reward %in% c("min","max")`
- `length(sig_th) == 1L && is.numeric`
- `length(window) == 1L && is.numeric && window > 0`
- `reset_on %in% c("best","any")`

---

## 5. 檔案讀取 — `read_sumstat()`

### 5.1 三種 `format` 預設

| `format` | `chrom` | `pos` | `id` | `value` | `test_filter` |
|---|---|---|---|---|---|
| `"plink2"` | `#CHROM` | `POS` | `ID` | `P` | `TRUE` |
| `"gpcm"` | `#CHROM` | `POS` | `ID` | `P_HPI` | `FALSE` |
| `"custom"` | (使用者必填) | (必填) | (必填) | (必填) | `FALSE` |

使用者顯式傳入的參數會 override 預設 (`%||%`)。

### 5.2 處理流程

1. **`match.arg(format)`** — 允許縮寫 (例如 `"p"` 即 `"plink2"`)
2. **`data.table::fread()`** 讀檔，自動偵測分隔符，支援壓縮；以 `tryCatch` 包裝錯誤訊息
3. **欄名正規化**：`#CHROM` → `CHROM` (避免 R 對 `#` 開頭欄名的引號麻煩)；`chrom_col` 同步更新
4. **欄位存在性檢查**：缺欄即 `stop()`
5. **TEST 過濾**：若 `test_filter = TRUE`
   - 若 `test_col` 不存在 → `warning()` 並跳過 (不終止)
   - 否則保留 `as.character(df[[test_col]]) == test_val` 之列，輸出 `message()` 顯示前後筆數
   - 過濾後 0 筆即 `stop()`
6. **染色體過濾**：若 `chrom` 非 `NULL`，以 `as.character()` 比對保留
7. **附加介面欄位**：`df$position`、`df$value` 經 `suppressWarnings(as.numeric())` 轉換
8. **NA 丟棄**：刪除 `position` 或 `value` 為 NA 之列，並 `message()` 報告丟棄筆數
9. **依 `position` 升冪排序** (`physical_merge()` 假設輸入已排序)
10. **建議 `reward`**：當 `value_col ∈ {LOG10_P, T_STAT, Z_STAT, CHISQ, F_STAT, T_STAT_Direct, T_STAT_TE, HPI}` → `"max"`，否則 `"min"`；`LOG10_P` 額外發出 `message()`

### 5.3 回傳結構

```r
list(
  data   = <data.frame，含原始欄位 + position + value，已排序>,
  reward = "min" | "max"
)
```

`reward` 的存在使得 `physical_merge(out$data, ..., reward = out$reward)` 不需要使用者再判斷方向。

---

## 6. 結果註解 — `annotate_blocks()`

### 6.1 用途
把 `physical_merge()` 簡化後的 block 表，與原始輸入 (含 BETA/SE/REF/ALT/...) 做 left join，每個 block 取代表性 SNP 那一列的完整欄位。

### 6.2 流程

1. **早退**：`nrow(blocks) == 0` 直接 return
2. **染色體欄解析**：使用者 `chrom_col` → `"CHROM"` → `"#CHROM"` → `NULL`
   (容忍未經 `read_sumstat()` 處理、保留 `#CHROM` 的輸入)
3. **ID 欄解析**：使用者 `id_col` → `"ID"` → `"SNP"` → `NA`
4. **主路徑 — 直接用列號**：若 `attr(blocks, "rps_row")` 存在、長度符合、且 `data$position[rps_row]` 與 `blocks$rps_BP` 完全相符，就以 `data[rps_row, ]` 取出代表列，用 `cbind` 接上 blocks。這是唯一能在同 bp 多 SNP (multi-allelic) 時標對 lead SNP 的方式。
5. **回退路徑 — 座標比對** (blocks 非本版 `physical_merge()` 產生時)：
   - 去重：有 chrom `!duplicated(data[, c(chrom_col, "position")])`；無 chrom `!duplicated(data$position)`
   - 代表列過濾：`paste(CHROM, rps_BP, sep=":")` 複合鍵或 `position %in% rps_BP`
   - **此路徑的限制**：同 bp 多 SNP 只保留首次，可能標到不顯著的那顆
6. 建立 `rps_BP = position`、`rps_ID = data[[id_col]]` (若可得)
7. 挑選欲攜帶的欄位：`c(rps_BP, [CHROM], [rps_ID], <其他原始欄位>)`，剔除 `position/value/id_col` 避免重複
8. **Join**：`merge(blocks, repr, by = merge_keys, all.x = TRUE)`，`merge_keys = c("CHROM","rps_BP")` 或 `"rps_BP"`
9. **欄位順序控制** by `keep_*` flags：
   - `meta = [serial?, CHROM?, start?, end?, rps_BP?, rps_ID?, rps_value?]`
   - `rest = setdiff(names(out), meta ∪ {標準 meta 欄})`
   - 最終 `out = out[, c(meta, rest)]`，依 `serial` 排序

### 6.3 `keep_*` 設計理由
為 downstream 工具 (例如 PRS pipeline、UKB-RAP Swiss Army Knife) 提供最小欄位輸出。注意：若後續還要用 `export_snp_list()`，`keep_rps_BP` 應保持 `TRUE` (fallback 識別碼)。

---

## 7. 匯出 — `export_snp_list()`

### 7.1 兩種模式

| `by_chrom` | 輸出 | 條件 |
|---|---|---|
| `FALSE` (預設) | 單一 `.txt` (每行一個 ID) | 無 |
| `TRUE` | `.zip`，內含 `snp_ch<CHROM>.txt` 多檔 | 需要 `CHROM` 欄 |

### 7.2 ID 欄位 fallback
`id_col` 預設依序 `rps_ID` → `rps_BP`，使得即使沒做 `annotate_blocks()`，也能用 bp 位置匯出 (用於 PLINK `--extract bed3` 之類)。

### 7.3 ZIP 子流程 (per-chrom)

```
path    <- normalizePath(path, mustWork = FALSE)     # 在 setwd 前先解析相對路徑
tmp_dir <- tempfile("physmerge_export_")
dir.create(tmp_dir)
old_wd  <- getwd()
on.exit({ setwd(old_wd); unlink(tmp_dir, recursive = TRUE) }, add = FALSE)

for ch in sort(unique(CHROM)):
    writeLines(ids[CHROM == ch], paste0("snp_ch", ch, ".txt"))

setwd(tmp_dir)
utils::zip(path, files = list.files(tmp_dir), flags = "-j")  # -j: 不存路徑
```

關鍵設計：
- `normalizePath()` **先於** `setwd()` 解析，避免 zip 寫到錯目錄
- `on.exit()` 同時負責**還原 wd** 與**清理 tmp dir** (即使 zip 失敗也安全)
- `-j` flag 讓 zip 內檔案不含目錄結構

### 7.4 早退與例外
- `nrow(blocks) == 0` → `warning()` + 早 return
- `id_col` 不存在 → `stop()`
- `by_chrom = TRUE` 但無 `CHROM` 欄 → `stop()`

---

## 8. 顯著性語意 (Significance semantics)

兩個 helper 將 p-value/test-statistic 的雙向比較抽象化：

```r
.is_significant(val, sig_th, reward)      # 是否達顯著閾值
  reward == "min": val <  sig_th
  reward == "max": val >  sig_th

.is_more_significant(val, best, reward)   # 是否比 best 更顯著
  reward == "min": val <  best
  reward == "max": val >  best
```

此抽象讓核心迴圈只需呼叫這兩個函式，**單一程式碼路徑同時處理 p-value 與 test statistic**。

---

## 9. 複雜度 (Complexity)

| 階段 | 時間 | 空間 |
|---|---|---|
| Sort by position | O(n log n) | O(n) |
| Forward scan | O(n) | O(n) — 預先配置 `out_*` 向量 |
| Collapse | O(b) (b = raw block count) | O(b) |
| Trim | O(b) | in-place |
| Per-chrom 分流 | O(n) extra (subset) | O(n) extra (lapply 結果) |

`out_*` 預配置長度 `n` 是保守上界 (每筆 SNP 最多開一個 block)；最後用 `[seq_len(block_count)]` 切回實際長度。

---

## 10. 邊界案例 (Edge cases) 與其處置

| 情境 | 處置位置 | 行為 |
|---|---|---|
| `nrow(data) == 0` | `.physical_merge_single` | 回傳 0-row template data frame (含正確欄名/型別) |
| 整段無顯著 SNP | scan 結束時 `block_count == 0` | 同上，回傳空表 |
| 結尾在 block 內 | 迴圈外 `if (in_block) close_block(last_pos)` | 正常收尾 |
| 單一 block 結果 | `.collapse_blocks` 早退 | 直接回傳，跳過 collapse |
| 跨染色體無欄位且 range > 250 Mb | `physical_merge()` top-level | `warning()` 但繼續執行 |
| `test_col` 不存在但 `test_filter = TRUE` | `read_sumstat` | `warning()` 並跳過，不 stop |
| `LOG10_P` 作為 value | `read_sumstat` | `message()` 並建議 `reward = "max"` |
| 同 bp 多 SNP | `annotate_blocks` | 只保留首次出現 (已知限制) |

---

## 11. 已知限制與未來方向 (Known Limitations / Future Work)

1. **同 bp 多 SNP**：已修正。`physical_merge()` 記錄代表 SNP 的輸入列號於 `attr(blocks, "rps_row")`，`annotate_blocks()` 據此取列。僅在 blocks 不帶該 attribute 時才退回舊的去重比對。實測某染色體 22 的 PLINK2 檔，542 個 block 中原有 2 個標錯。
2. **`reset_on = "any"` 的代表 SNP 行為**：window 重置與代表更新分離，是刻意設計。若使用者期望「最後出現的顯著 SNP」當代表，需另外擴充。
3. **無 LD 資訊**：物理距離合併在高 LD 區 (如 MHC) 可能過度或不足合併。建議搭配 region-specific window。
4. **跨染色體保護**：250 Mb 啟發式僅在沒有 `CHROM` 欄時觸發。最佳實踐：永遠提供 `CHROM`。
5. **無 header 檔案**：`fread(header = TRUE)` 寫死，無 header 之檔案需先預處理。
6. **`data.table` 為 Suggests**：若使用者未安裝，`read_sumstat()` 會失敗。可考慮改成 hard Import 或加 fallback `utils::read.table`。
7. **`shiny`/`DT`/`shinyjs`** 列於 Suggests 但目前未被任何 exported 函式使用，可能是預留 UI。

---

## 12. 測試 (Tests)

`tests/testthat/` 目錄存在 (具體案例未在本 spec 範圍內列舉)。建議覆蓋：
- 三 passes 的 invariant (no overlap、serial 連續性)
- `reset_on` 兩模式在同一輸入下的差異
- per-chrom 分流的 serial 重編
- `read_sumstat` 三種 format 的欄位 fallback
- `export_snp_list` 的 wd 復原 (on.exit) 行為
- 同 bp 兩顆 SNP 時代表 SNP 的正確性，以及無 `rps_row` attribute 時的回退
- 邊界：空資料、單一 block、跨 chr 邊界

---

## 13. 函式呼叫流程圖 (End-to-end pipeline)

```
            ┌────────────────────────┐
   file ──▶ │ read_sumstat()         │
            │  - fread               │
            │  - normalize #CHROM    │
            │  - TEST filter         │
            │  - chrom filter        │
            │  - append position/value
            │  - drop NA, sort       │
            │  - suggest reward      │
            └──────────┬─────────────┘
                       │ list(data, reward)
                       ▼
            ┌────────────────────────┐
            │ physical_merge()       │
            │  - validate            │
            │  - per-chrom dispatch  │
            │  → .physical_merge_single
            │     · Pass 1: forward scan
            │     · Pass 2: collapse │
            │     · Pass 3: trim     │
            └──────────┬─────────────┘
                       │ blocks (serial, CHROM, start, end, rps_BP, rps_value)
                       ▼
            ┌────────────────────────┐
            │ annotate_blocks()      │
            │  - dedup on (CHROM, position)
            │  - filter to rep rows  │
            │  - merge back columns  │
            │  - apply keep_* flags  │
            └──────────┬─────────────┘
                       │ annotated blocks (+ BETA/SE/REF/ALT/...)
                       ▼
            ┌────────────────────────┐
            │ export_snp_list()      │
            │  - by_chrom = FALSE → .txt
            │  - by_chrom = TRUE  → .zip (per chrom)
            └──────────┬─────────────┘
                       ▼
                  下游分析 (PRS、LDSC、coloc、UKB-RAP)
```

---

## 附錄 A — 主要符號 (Glossary)

| 符號 | 意義 |
|---|---|
| `sig_th` | 顯著閾值 (e.g. 5e-8) |
| `window` | 視窗大小 (bp，e.g. 500000) |
| `steps` | 當前 block 內視窗剩餘 bp 餘額 |
| `reward` | `"min"` (p-value) 或 `"max"` (test stat) |
| `reset_on` | `"best"` 或 `"any"` — 視窗重置觸發條件 |
| `rps_BP` / `rps_value` / `rps_ID` | Representative SNP 之位置 / 值 / ID |
| `serial` | Block 連續編號 (跨 chr 重編) |
