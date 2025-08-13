# Molecular Property Prediction — Neural Networks

This project explores methods to improve a **chemical language model** for predicting **lipophilicity** from molecular **SMILES** strings. It is organized into three tasks:

* **Task 1:** Supervised fine-tuning for regression and **unsupervised MLM pre-finetuning** (masked-token prediction) to boost downstream performance.&#x20;
* **Task 2:** **Influence function**–based data selection (with **LiSSA** inverse-HVP approximation) to choose the most impactful external samples for training.&#x20;
* **Task 3:** **Parameter-efficient fine-tuning (PEFT)** strategies — **BitFit**, **LoRA**, and **IA³** — to reduce compute while preserving accuracy.&#x20;

> **NNTI Project (WiSe 2025)** — Abdelsalam Helala , Mostafa Hammouda , Loulwah Arnaout .&#x20;

---

## Repository structure

```
Neural-Networks/
└── Molecular Property Prediction/
    ├── notebooks & scripts for Task 1–3
    └── report/ (LaTeX/PDF)
```

---

## Environment

* Python 3.10+
* PyTorch, Hugging Face Transformers
* scikit-learn, matplotlib
* (Task 3) PEFT libraries (e.g., LoRA/IA³ implementations)
* (Optional) Weights & Biases for logging

Install (example):

```bash
pip install torch transformers scikit-learn matplotlib wandb
# If using PEFT implementations:
pip install peft
```

---

## Data

* **Lipophilicity dataset**: SMILES + target lipophilicity values.
* Data are split **80% train / 20% test**; batch size **16** in baseline experiments.&#x20;

---

## Methods

### Task 1 — Baseline fine-tuning & MLM pre-finetuning

1. Load a **pretrained chemical LM** (e.g., MoLFormer) and add a **regression head** (MSE loss, Adam/AdamW).
2. Train for up to **20 epochs**, early-stopping around **8 epochs** when eval loss stabilizes (overfitting beyond that).&#x20;
3. **Unsupervised pre-finetuning** via **MLM** on SMILES (cross-entropy) improves downstream smoothness of the loss curve and test metrics.&#x20;

### Task 2 — Influence function data selection

* Compute per-sample influence using **LiSSA** to approximate the **inverse Hessian–vector product**; rank external samples by influence and **select top 50%** for further training. This **improves generalization** vs. Task-1 baselines.&#x20;

### Task 3 — Parameter-efficient fine-tuning (PEFT)

* **BitFit:** update **biases only** (fastest, lowest memory).
* **LoRA:** low-rank adapters in attention (here **r = 8**, **α = 16**).
* **IA³:** multiplicative scaling of activations (lightweight).
* Example training schedule:

  * BitFit: **30** epochs, **2e-5** LR
  * LoRA: **20** epochs, **1e-5** LR
  * IA³: **20** epochs, **5e-6** LR
    BitFit was most efficient and achieved the best eval loss in our runs; LoRA offered flexibility at higher runtime; IA³ was light but less impactful.&#x20;

---

## Results (summary)

| Model / Strategy                            | MSE        | MAE        | R²         |
| ------------------------------------------- | ---------- | ---------- | ---------- |
| Baseline fine-tune                          | 0.4723     | 0.5316     | 0.6804     |
| + **Unsupervised MLM pre-finetuning**       | 0.4259     | 0.4942     | 0.7117     |
| **Influence-selected data (IFDS)** (Task 2) | **0.3945** | **0.4772** | **0.7330** |

Unsupervised MLM pre-finetuning and IF-based selection both **lower error** and **raise R²** over the baseline; IFDS performs best overall.&#x20;

---

## Reproducibility

* Fix random seeds; log training & evaluation curves (e.g., WandB).
* Monitor **train vs. eval loss** to avoid overfitting (Task 1 stabilized \~epoch 8).&#x20;
* For influence functions, ensure stable LiSSA settings (damping/scale, iterations) before ranking samples.&#x20;

---

## References (selection)

Influence functions (Koh & Liang), LiSSA (Agarwal et al.), Transformers (Wolf et al.), Adam (Kingma & Ba), matplotlib (Hunter), scikit-learn (Pedregosa). Full citations are provided in the report.&#x20;
