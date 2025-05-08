import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# 定义颜色名称与十六进制码
colors = [
    ("Deep Blue", "#183F5E"),
    ("Light Blue", "#21639A"),
    ("Yellow", "#F6D565"),
    ("Green", "#3EAEA4"),
    ("Red", "#EC553C"),
]

fig, ax = plt.subplots(figsize=(7, 2))
for idx, (name, hexcode) in enumerate(colors):
    rect = Rectangle((idx, 0), 1, 1, color=hexcode)
    ax.add_patch(rect)
    ax.text(idx + 0.5, -0.1, f"{name}\n{hexcode}", 
            ha='center', va='top', fontsize=10)

ax.set_xlim(0, len(colors))
ax.set_ylim(-0.2, 1)
ax.axis('off')
plt.tight_layout()
plt.show()