import re
import matplotlib.pyplot as plt

def parse_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    thresholds = [float(t) if 'e' in t else float(t)
                  for t in re.findall(r"Threshold: ([\de\.-]+)", content)]
    mcc_values = [float(m) for m in re.findall(r"MCC: ([\d\.]+)", content)]

    # Filtro: threshold <= 1e-5
    filtered = [(t, m) for t, m in zip(thresholds, mcc_values) if t <= 1e-5]
    return zip(*filtered)  # restituisce (thresholds, mcc_values)

# Percorsi ai file
set1_path = 'performance_set1_thresholds.txt'
set2_path = 'performance_set2_thresholds.txt'

# Parsing dei dati filtrati
thresholds1, mcc1 = parse_file(set1_path)
thresholds2, mcc2 = parse_file(set2_path)

# Plot
plt.figure(figsize=(10, 6))
plt.plot(thresholds1, mcc1, marker='o', color='blue', label='Set 1')
plt.plot(thresholds2, mcc2, marker='s', color='red', label='Set 2')
plt.xscale('log')
plt.ylim(0.98, 1.001)
plt.xlabel('Threshold')
plt.ylabel('MCC')
plt.title('Confronto MCC tra Set 1 e Set 2 (threshold ≤ 1e-5)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
