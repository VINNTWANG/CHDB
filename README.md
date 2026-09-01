# CHD Gene Family Variant & HPO Explorer
# CHD 基因家族变异与 HPO 表型探索器

[English](#english) | [中文](#chinese)

---

<a name="english"></a>
## English Version

This project is a Streamlit-based interactive dashboard designed to explore and analyze ClinVar variants and Human Phenotype Ontology (HPO) annotations for the CHD (Chromodomain Helicase DNA-binding) gene family (CHD1-CHD9).

### 🚀 Features

#### 1. Gene Explorer (ClinVar)
- **Positional Plots:** Interactive Lollipop, Waterfall, and Density plots to visualize variant distribution.
- **Variant Analysis:** Breakdowns of variant types vs. clinical significance.
- **Phenotype Analysis:** Relationships between variant types and ClinVar reported phenotypes.
- **Origin Analysis:** Distribution of variant origins (e.g., de novo, inherited).

#### 2. CHD Family Comparison
- **Visualizations:** Pathogenic variant heatmaps, significance distribution across the family, and Sankey diagrams linking types to phenotypes.
- **Summary Statistics:** Protein lengths, total variant counts, and pathogenic variant ratios for CHD1-9.

#### 3. HPO Phenotype Explorer
- **Overview:** Frequency of HPO terms and associated disease links.
- **Visuals:** Phenotype category breakdown, Word Clouds, and Gene-Phenotype-Disease network graphs.

### 📊 Data Sources
- **ClinVar:** Variant data (`CHDx_clinvar.txt`).
- **UniProt:** Protein domain info (`CHDx_domains.json`).
- **HPO:** Phenotype annotations (`chd_phenotype_data/CHDx_annotations.json`).

### 🛠️ Setup and Running
1. **Install Dependencies:**
   ```bash
   pip install -r requirements.txt
   ```
2. **Run Application:**
   ```bash
   streamlit run app11.py
   ```

---

<a name="chinese"></a>
## 中文版

本项目是一个基于 Streamlit 的交互式仪表板，旨在探索和分析 CHD（染色质解旋酶 DNA 结合）基因家族（CHD1-CHD9）的 ClinVar 变异及人类表型本体（HPO）注释数据。

### 🚀 功能特性

#### 1. 基因浏览器 (ClinVar)
- **位置图谱：** 交互式棒棒糖图 (Lollipop)、瀑布图 (Waterfall) 和密度图，可视化变异在蛋白质序列上的分布。
- **变异分析：** 变异类型与临床意义（良性/致病等）的详细分类。
- **表型分析：** 变异类型与 ClinVar 报告表型之间的关联。
- **来源分析：** 变异来源分布（如：新发变异 de novo、遗传性变异）。

#### 2. CHD 家族对比
- **可视化对比：** 致病性变异热图、家族成员变异意义分布对比、以及连接变异类型与表型的桑基图 (Sankey)。
- **统计摘要：** CHD1-9 的蛋白质长度、变异总数及致病变异比例概览。

#### 3. HPO 表型探索器
- **表型概览：** HPO 术语频率及关联疾病统计。
- **视觉呈现：** 表型类别划分、词云图 (Word Cloud) 以及“基因-表型-疾病”交互网络图。

### 📊 数据来源
- **ClinVar:** 变异数据 (`CHDx_clinvar.txt`)。
- **UniProt:** 蛋白质结构域信息 (`CHDx_domains.json`)。
- **HPO:** 表型注释数据 (`chd_phenotype_data/CHDx_annotations.json`)。

### 🛠️ 安装与运行
1. **安装依赖：**
   ```bash
   pip install -r requirements.txt
   ```
2. **启动应用：**
   ```bash
   streamlit run app11.py
   ```

---

### 👨‍💻 Developer / 开发人员
**Soap**

**Contact:** [soap628@proton.me](mailto:soap628@proton.me)
