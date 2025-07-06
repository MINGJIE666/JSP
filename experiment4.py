import matplotlib.pyplot as plt
import pandas as pd

# 设置Palatino字体
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

# 读取数据
file='D:\\工作备份\\论文\\JSPJOC\\table\\experiment4.xlsx'
d = pd.read_excel(file)

x = d["a"]
y1 = d["b"]
y2 = d["c"]
y3 = d["d"]
y4 = d["e"]

plt.figure(figsize=(20, 13)) #画布大小
plt.plot(x, y2, marker='^', ms=25,  lw=3, color='b', label='HTSFP4')
plt.plot(x, y3, marker='s', ms=25,  lw=3, color='g', label='HTSFP5')
#plt.plot(x, y4, marker='X', ms=25,  lw=3, color='y', label='HTSFP1$_{15}$')
plt.plot(x, y1, marker='o', ms=25, lw=3, color='r', label='HTSFP')


plt.ylim(-0.03, 1)
plt.xlim(-0, 12)
plt.xticks(rotation=-45) # 横坐标标签字体大小

# 设置坐标轴数值字体大小
plt.tick_params(axis='both', which='major', labelsize=35)


plt.grid(True, axis='x', linestyle='--', lw=0.5) 
plt.legend(fontsize=35, loc=2, edgecolor='black') 
plt.xlabel('Instance', fontsize=50, fontweight='bold')
plt.ylabel('Gap to the Best-known Solutions (\%)', fontsize=50)
plt.gca().set_axisbelow(True)

plt.tight_layout()
plt.savefig('D:\\工作备份\\论文\\JSPJOC\\table\\experiment4.eps')
plt.savefig('D:\\工作备份\\论文\\JSPJOC\\table\\experiment4.pdf')
plt.show()