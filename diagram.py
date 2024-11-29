import matplotlib.pyplot as plt
import numpy as np

# Initialize the figure and axes
fig, ax = plt.subplots(2, 1, figsize=(10, 12))

# Upper plot: multi-sheeted complex plane L = C
ax[0].set_xlim(-2, 4)
ax[0].set_ylim(-10, 10)
ax[0].set_aspect("equal")
ax[0].set_title("L = ℂ (Covering space)", fontsize=14)

# Grid points (ln(2) + k*2πi)
ln2 = np.log(2)
k_values = np.arange(-2, 3)
for k in k_values:
    ax[0].scatter([ln2], [2 * np.pi * k], color="black", zorder=5)
    ax[0].text(ln2 + 0.2, 2 * np.pi * k, f"ln2+i{k}2π", fontsize=10)

# Vertical lines
ax[0].plot([ln2, ln2], [-4 * np.pi, 4 * np.pi], color="purple", lw=2, zorder=3, label="Periodic jumps")
ax[0].plot([0, 0], [-4 * np.pi, 4 * np.pi], color="red", lw=2, linestyle="--", label="Reference line")

# Horizontal plane shading for 0th sheet
ax[0].fill_betweenx([-np.pi, np.pi], -2, 4, color="pink", alpha=0.3, label="0th Sheet")

# Annotations and labels
ax[0].annotate("", xy=[ln2, 2 * np.pi], xytext=[ln2, 0], arrowprops=dict(arrowstyle="->", color="blue"))
ax[0].annotate("", xy=[ln2, -2 * np.pi], xytext=[ln2, 0], arrowprops=dict(arrowstyle="->", color="red"))
ax[0].text(ln2, -5, "Multi-sheeted", fontsize=10)
ax[0].legend()

# Middle mapping label
ax[0].annotate(r"$\eta(z)=e^z-1$", xy=(1, -7), xytext=(2, -8), fontsize=14, arrowprops=dict(arrowstyle="->"))

# Lower plot: punctured complex plane M = ℂ \ {-1}
ax[1].set_xlim(-3, 3)
ax[1].set_ylim(-3, 3)
ax[1].set_aspect("equal")
ax[1].set_title(r"M = ℂ \ {-1} (Punctured plane)", fontsize=14)

# Draw puncture point at -1
ax[1].scatter([-1], [0], color="black", s=100)
ax[1].text(-1.2, 0.2, "-1", fontsize=12)

# Circles representing preimages under eta(z)
theta = np.linspace(0, 2 * np.pi, 100)
for r, color in zip([0.5, 1, 1.5], ["red", "blue", "orange"]):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    ax[1].plot(x, y, color=color, label=f"r = {r}")

# Annotations
ax[1].annotate("", xy=[-1, 0], xytext=[-1, -2], arrowprops=dict(arrowstyle="->", color="purple"))
ax[1].text(-1.5, -1, "Mapping loops", fontsize=10)
ax[1].legend()

# Display the plot
plt.tight_layout()
plt.show()
