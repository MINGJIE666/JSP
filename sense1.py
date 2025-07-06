import matplotlib.pyplot as plt
import pandas as pd

# 设置Palatino字体
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

# 读取数据
file='sense1.xlsx'
d = pd.read_excel(file)

x = d["a"]
y1 = d["b"]
y2 = d["c"]


plt.figure(figsize=(10, 10)) #画布大小
plt.plot(x, y1, marker='s', ms=20,  lw=3, color='r', label='$\Phi_{best}$')
plt.plot(x, y2, marker='^', ms=25,  lw=3, color='b', label='$\Phi_{avg}$')



plt.ylim(-0.1, 1.1)
plt.xlim(900, 5100)
#plt.xticks(rotation=-45) # 横坐标标签字体大小

# 设置坐标轴数值字体大小
plt.tick_params(axis='both', which='major', labelsize=35)


plt.grid(True, axis='x', linestyle='--', lw=0.5) 
plt.legend(fontsize=35, loc=2, edgecolor='black') 
plt.xlabel('Parameter Valve', fontsize=35, fontweight='bold')
plt.ylabel('Gap(\%)', fontsize=35)
plt.gca().set_axisbelow(True)

plt.tight_layout()
plt.savefig('sense1.eps')
plt.savefig('sense1.pdf')
plt.show()