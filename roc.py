import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

# File di input aggiornati
files = ["set_1.class", "set_2.class"]
y_true = []
scores = []

# Leggi entrambi i file
for fname in files:
    with open(fname) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 4:
                try:
                    label = int(parts[1])       # Etichetta dalla seconda colonna
                    score = float(parts[2])     # E-value dalla terza colonna
                    scores.append(-score)       # Più piccolo l’E-value = più confidenza
                    y_true.append(label)
                except ValueError:
                    continue  # Salta righe malformate

# Calcolo ROC
fpr, tpr, thresholds = roc_curve(y_true, scores)
roc_auc = auc(fpr, tpr)

# Plot della curva ROC
plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, label=f"ROC curve (AUC = {roc_auc:.2f})", color='darkgreen')
plt.plot([0, 1], [0, 1], linestyle='--', color='gray')
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve for HMM Classifier")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("roc_curve.png")
plt.show()
