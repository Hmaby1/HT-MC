import numpy as np
import matplotlib.pyplot as plt
import math

# Boltzmann constant
k = 8.617333262145E-5
T = 200

# 创建delta_E
x = np.linspace(0,1, 100)

# 计算 y 数据，y = exp(-x / (cls.k * T))
y = np.exp(-x / (k * T))

# 绘制函数图像
plt.plot(x, y,         #折线图
        color="steelblue",
        linewidth=1,
        )

# 填充函数与 x 轴之间的区域（阴影）
plt.fill_between(x, y, 
                 color='lightblue', 
                 alpha=0.4,)
# 2. xlabel, ylabel
plt.xlabel("delta_E", fontsize=16, fontweight="bold")
plt.ylabel("Acceptance", fontsize=16, fontweight="bold")

# 3. xaxis, yaxis
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

legend_font = {"size" : 13, 
        "weight": "bold"
        }
plt.legend(loc=(0.6, 0.75),
        prop=legend_font,
        frameon=False)

plt.show()
