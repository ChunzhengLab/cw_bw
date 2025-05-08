import pandas as pd
import matplotlib.pyplot as plt

# 1. 读取 CSV 文件
df = pd.read_csv('scan_betascale.csv', sep=',')

# 2. 确保数值列类型正确
df['betaT']     = pd.to_numeric(df['betaT'],     errors='coerce')
df['v2']        = pd.to_numeric(df['v2'],        errors='coerce')
df['v2_Err']    = pd.to_numeric(df['v2_Err'],    errors='coerce')
df['centrality'] = df['centrality'].astype(str)  # 作为分类标签

# 3. 按 centrality 分组绘图（带误差条）
plt.figure(figsize=(8,6))
for cent, group in df.groupby('centrality'):
    group_sorted = group.sort_values('betaT')
    plt.errorbar(
        group_sorted['betaT'],
        group_sorted['v2'],
        yerr=group_sorted['v2_Err'],
        marker='o',
        linestyle='-',
        capsize=3,           # 误差条横线长度
        label=f'Centrality {cent}'
    )

# 4. 美化
plt.xlabel(r'$\beta_T$')
plt.ylabel(r'$v_2$')
plt.title(r'$v_2$ vs. $\beta_T$ for different centralities')
plt.legend(title='Centrality')
plt.grid(True)

# 5. 保存与显示
plt.tight_layout()
plt.savefig('v2_vs_betaT_with_errors.pdf')
plt.show()