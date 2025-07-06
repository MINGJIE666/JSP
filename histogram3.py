import matplotlib.pyplot as plt
import numpy as np

# 设置字体为 TeX 的默认字体
plt.rc('font', family='serif')
plt.rc('text', usetex=True)  # 启用 LaTeX 渲染

# 数据
categories = ['High', 'Random', 'Low']
adjacency = [57, 21, 18]  # adjacency的比例
precedence = [89, 31, 26]  # precedence的比例

x = np.arange(len(categories))  # 横坐标位置
bar_width = 0.25  # 柱子宽度

# 绘图
fig, ax = plt.subplots(figsize=(4, 3))

# 绘制柱状图
bars1 = ax.bar(x , adjacency, bar_width, label='Adjacent Job Order', color='grey')
#bars2 = ax.bar(x + bar_width / 2, precedence, bar_width, label='Precedence Order', color='blue')

# 添加标签和标题
ax.set_xlabel('Solution Quality', fontsize=12)
ax.set_ylabel('Proportion(\%)', fontsize=12)
ax.set_xticks(x)
ax.set_xticklabels(categories)
ax.set_ylim(0, 100)

# 添加图例
#ax.legend()

# 添加柱子上方的数值标注
for bars in [bars1]:
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, height + 2, f'{height}%', ha='center', va='bottom', fontsize=10)

# 显示图表
plt.tight_layout()
plt.savefig('histogram3.pdf', format='pdf')
plt.show()
