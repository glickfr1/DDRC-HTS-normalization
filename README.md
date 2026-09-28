# DDRC High-Throughput Screening Normalization Tool

Developed by Chloe Larson, Francesca Curreli, Loreto Carvallo-Torres, and J. Fraser Glickman for drug-discovery and high-throughput screening workflows at the Fisher Drug Discovery Resource Center (DDRC), The Rockefeller University.

Copyright © 2026 The Rockefeller University.

## License

This software is licensed under the **Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License (CC BY-NC-ND 4.0)**.

For licensing inquiries or requests for permissions beyond the terms of this license, please contact The Fisher Drug Discovery Resource Center by opening an Issue in this GitHub repository.
## What does this tool do?

This tool takes raw high-throughput screening (HTS) measurements and a separate compound database, links the two files by **plate + well**, and produces normalized screening results.

It calculates:

- **Percent inhibition** for each well
- **Raw data z-score** for each well (plate-based)
- **Percent inhibition z-score** for each well (plate-based)
- **Run percent inhibition z-score** for each well (run-based)
- **Z′ factor** for each plate
- plate QC classification
- potential hits above a selected percent-inhibition threshold

It produces CSV, Excel, PDF, and Parquet output. Parquet output requires `pyarrow`.

## The important idea: two files work together

An HTS instrument usually reports **what happened in a well**, but the instrument data may not contain the chemical identity of the compound in that well.

For example, the screening data might contain:

| Plate | Well | Fluorescence (RFU) |
|---|---|---:|
| 1827 | G16 | 6250 |

A separate compound database might contain:

| Plate | Well | Molecule Name | Structure (CXSMILES) |
|---|---|---|---|
| 1827 | G16 | RU-0000353 | `COC1=C(OC)C(OC)=CC(CC2=C(N)N=C(N)N=C2)=C1` |

The program connects these records using **Plate + Well**.

This separation is useful because the screening measurements and compound information can be maintained independently while still being joined into one analysis-ready dataset.

## What is a SMILES or CXSMILES string?

SMILES is a text representation of a chemical structure. CXSMILES is an extended form of SMILES that can include additional chemical information.

The normalization program does **not** calculate molecular descriptors or otherwise analyze the structures. It carries the selected compound identifier and structure field from the compound database into the normalized output.

This means the resulting file can be opened in chemistry software such as DataWarrior, which can read structures and calculate molecular properties and descriptors from them.

## Input files

### Screening dataset CSV

The screening dataset must contain:

- a **plate identifier** column
- a **well identifier** column
- the **raw data/readout** column you want to analyze

The program asks you to select these columns by number, so the columns do not have to have particular names.

For example:

```text
Plate,Well,Control State,Fluorescence (RFU)
1827,G16,,6250
1818,G14,,4040
2033,G04,,5250
```

An existing control-state column may be present in the input data, but the control layout selected when the program is run determines the standardized `Control State` used for normalization.

### Compound database CSV

The compound database must contain:

- a **plate identifier** column
- a **well identifier** column
- a **compound identifier** column
- a **structure** column containing SMILES, CXSMILES, or another structure representation that you want carried into the output

The program asks you to select these columns by number.

Optional molecule/name and batch columns can also be selected.

For example:

```text
Molecule Name,Structure (CXSMILES),Batch Name,Plate,Well
RU-0000353,COC1=C(OC)C(OC)=CC(CC2=C(N)N=C(N)N=C2)=C1,AA-008,1827,G16
```

## Control wells

The program supports three common control layouts:

1. **Entire column(s)** — for example, columns 23 and 24
2. **Entire row(s)** — for example, rows A and B
3. **Specific wells** — for example, A01,A02,A03

For each layout, you specify which locations contain the **POSITIVE** control and which contain the **NEGATIVE** control.

The program adds a standardized `Control State` column to the output with one of three values:

- `POSITIVE`
- `NEGATIVE`
- blank for experimental wells

## Percent inhibition

Percent inhibition is calculated separately for each plate using the mean signal of the negative and positive control populations:

```text
Percent inhibition =
100 - ((well value - positive-control mean)
       / (negative-control mean - positive-control mean) × 100)
```

The interpretation depends on the assay design and the assignment of positive and negative controls. In the usual inhibition format, values near 0% resemble the negative-control signal and values near 100% resemble the positive-control signal.

## Understanding the z-scores

The program reports three different z-score measurements. They answer different questions because they use different reference populations.

### Raw data z-score — plate-based

The **Raw data z-score** standardizes the selected raw readout within each individual plate:

```text
(raw value - plate mean) / plate standard deviation
```

It tells you how unusual a well's raw measurement is relative to the other measurements on that plate.

### Percent inhibition z-score — plate-based

The **Percent inhibition z-score** standardizes percent inhibition within each individual plate.

It tells you how unusual a well's normalized inhibition value is relative to the other wells on that plate.

### Run percent inhibition z-score — run-based

The **Run percent inhibition z-score** standardizes percent inhibition relative to the experimental wells across the entire screening run.

It tells you how unusual a well's normalized inhibition value is relative to the overall distribution of experimental wells in the run.

### Why are there three z-scores?

The distinction is the reference population:

| Output | Reference population | Question it answers |
|---|---|---|
| **Raw data z-score** | One plate | How unusual is the raw readout on this plate? |
| **Percent inhibition z-score** | One plate | How unusual is the normalized inhibition on this plate? |
| **Run percent inhibition z-score** | Entire run | How unusual is the normalized inhibition across the run? |

These are **well-level measurements**. They should not be confused with Z′, which is a plate-level assay-quality statistic.

## Z′ factor

Z′ is calculated separately for each plate from the negative- and positive-control populations:

```text
Z′ = 1 - 3 × (SDnegative + SDpositive)
         / |Meannegative - Meanpositive|
```

Current QC categories are:

- **GOOD:** Z′ ≥ 0.5
- **BORDERLINE:** 0 ≤ Z′ < 0.5
- **BAD:** Z′ < 0

Z′ is a **plate-level assay-quality measurement**. It is not a z-score and is not calculated independently for every well.

## Potential hits

A well is identified as a potential hit when it is an experimental well (not a control) and its percent inhibition is greater than or equal to the selected threshold.

The default threshold is **50% inhibition** and can be changed when the program is run.

## Running the program

### 1. Install Python

Install a current version of Python 3.

### 2. Install the required packages

Open a terminal in the tool directory and run:

```bash
pip install -r requirements.txt
```

The requirements include `pyarrow`, which is used to create Parquet output.

### 3. Run interactively

The easiest way to start is:

```bash
python hts_normalization.py
```

The program will ask you to select the screening and compound database files, map the relevant columns, define the control layout, choose the hit threshold, and specify the output directory.

Press Enter to accept a displayed default where one is provided.

### 4. Run with command-line options

Once you are comfortable with the tool, the same analysis can be run without prompts. For example:

```bash
python hts_normalization.py \
  --data "testdata.csv" \
  --database "testdb.csv" \
  --data-plate "Plate" \
  --data-well "Well" \
  --assay "Fluorescence (RFU)" \
  --db-plate "Plate" \
  --db-well "Well" \
  --db-compound "Molecule Name" \
  --db-structure "Structure (CXSMILES)" \
  --control-layout columns \
  --positive-spec 24 \
  --negative-spec 23 \
  --hit-threshold 50 \
  --output Output
```

## Example files

The `examples` directory contains screening data and a compound database that can be used to test the program.

```text
examples/
├── testdata.csv
└── testdb.csv
```

The files illustrate the basic relationship between screening data and compound data using **Plate + Well** as the join key.

## Output

The analysis produces:

- `*_normalized_*.csv` — full normalized dataset
- `*_normalized_*.parquet` — full normalized dataset in Parquet format
- `*_top_hits_*.xlsx` — potential hits and plate QC
- `*_plate_QC_*.csv` — plate-level Z′ and QC classification
- `*_plate_report_*.pdf` — plate heatmaps and QC information

The full normalized dataset retains the compound identifier and structure field selected from the compound database, along with the raw readout, normalization results, z-scores, hit status, and control state.

## Typical workflow

A simple HTS workflow using this tool is:

```text
Screening instrument data
          +
Compound database
          ↓
       Plate + Well
          ↓
       Normalization
          ↓
   Percent inhibition
          ↓
     z-scores + Z′
          ↓
      Hit selection
          ↓
 Analysis of compound structures
```

The normalized CSV or Parquet output can then be used for further analysis. Because the compound structure is retained as SMILES/CXSMILES, chemistry software such as DataWarrior can be used to visualize structures and calculate molecular descriptors without requiring those descriptors to be stored in the screening database.

## Current assumptions

The current version is designed around the DDRC screening workflow and assumes:

- 384-well plates (A–P × 1–24)
- plate + well identify a screening record
- positive and negative controls can be described as columns, rows, or specific wells
- screening and compound database records can be joined by plate + well

The input column names themselves are flexible because the program asks the user to select them.

If your screening workflow is different, the Python source can be modified.

## Need to modify it?

The program is intentionally provided as open Python source rather than as a black-box application.

If you need a modification—for example:

- a different plate format
- a different control arrangement
- a different input-file format
- another hit criterion
- additional calculations
- different report formats

you can ask an AI coding assistant such as ChatGPT to modify the script. A useful approach is to upload the Python file and describe the change you need.

**Always review and validate modifications before using the results for scientific decision-making.**

## DDRC

Developed for drug-discovery and high-throughput screening workflows at the Fisher Drug Discovery Resource Center (DDRC), The Rockefeller University.

