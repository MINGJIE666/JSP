import matplotlib.pyplot as plt
import pandas as pd

# 设置Palatino字体
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

# 读取数据
file='experiment3.xlsx'
d = pd.read_excel(file)

x = d["a"]
y1 = d["b"]
y2 = d["c"]
y3 = d["d"]
y4 = d["e"]

plt.figure(figsize=(20, 13)) #画布大小
plt.plot(x, y2, marker='^', ms=2,  lw=15, color='b', label='HTSFP2')
plt.plot(x, y3, marker='s', ms=2,  lw=15, color='g', label='HTSFP3')
plt.plot(x, y4, marker='s', ms=2,  lw=25, color='y', label='HTSFP4')
#plt.plot(x, y4, marker='X', ms=25,  lw=3, color='y', label='HTSFP1$_{15}$')
plt.plot(x, y1, marker='o', ms=2, lw=15, color='r', label='HTSFP')


plt.ylim(0, 100)
plt.xlim(-0, 10)
#0plt.xticks(rotation=-45) # 横坐标标签字体大小

# 设置坐标轴数值字体大小
plt.tick_params(axis='both', which='major', labelsize=35)


plt.grid(True, axis='x', linestyle='--', lw=0.5) 
plt.legend(fontsize=35, loc=5, edgecolor='black') 
plt.xlabel('Iteration $(10^6)$', fontsize=50, fontweight='bold')
plt.ylabel('Error Rate(\%)', fontsize=50)
plt.gca().set_axisbelow(True)

plt.tight_layout()
plt.savefig('experiment3.eps')
plt.savefig('experiment3.pdf')
plt.show()