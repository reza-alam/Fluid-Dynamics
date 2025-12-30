import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib.gridspec as gridspec

# %% 1. USER CONFIGURATION
# ---------------------------------------------------------
Bo = 100.0              # Bond Number
kappa_max = 25.0        # Max range for the plot x-axis

# Slider Configuration for Kappa_0
k0_min = 0.0
k0_max = 12.5
k0_step = 0.1
k0_init = 0.0

num_points = 2000       # Resolution
# ---------------------------------------------------------

k_grid = np.linspace(0, kappa_max, num_points)

# %% 2. PHYSICS FUNCTIONS

def get_omega(k, Bo):
   """
   Calculates dispersion relation:
   Omega = sqrt( (k + k^3/Bo) * tanh(k) )
   """
   term_gravity = k
   term_capillary = (k**3) / Bo
   term_depth = np.tanh(k)
   w_sq = (term_gravity + term_capillary) * term_depth
   return np.sqrt(np.maximum(0, w_sq))

omega_base = get_omega(k_grid, Bo)

# %% 3. FIGURE AND PLOT SETUP

# Use Matplotlib's internal Computer Modern font for math
plt.rcParams.update({
   "text.usetex": False,
   "font.family": "serif",
   "mathtext.fontset": "cm",     
   "axes.unicode_minus": False
})

fig = plt.figure(figsize=(10, 7), facecolor='white')
fig.canvas.manager.set_window_title('Interactive Dispersion Relation')

gs = gridspec.GridSpec(2, 1, height_ratios=[1, 0.25], hspace=0.3)
ax = fig.add_subplot(gs[0])

# --- LABELING ---
ax.set_xlabel(r'$\kappa$', fontsize=16)
ax.set_ylabel(r'$\Omega$', fontsize=16)

# FIX: Use \mathbf instead of \textbf
title_str = r'$\mathbf{Dispersion\ Relation}\ (\mathrm{Bo} = %.1f)$' % Bo
ax.set_title(title_str, fontsize=16)

ax.grid(True, linestyle=':', alpha=0.6)
ax.set_xlim(0, kappa_max)
ax.set_ylim(0, np.max(omega_base))

# --- Initial Plot Objects ---
line_base, = ax.plot(k_grid, omega_base, 'b-', linewidth=2, label=r'$\Omega(\kappa)$')
line_shift, = ax.plot([], [], 'r--', linewidth=2, label=r'$\Omega(\kappa - \kappa_0) + \Omega_0$')
point_k0, = ax.plot([], [], 'ro', markersize=8, zorder=5)
point_int, = ax.plot([], [], 'go', markersize=10, zorder=6, label=r'Intersection')

ax.legend(loc='upper left', fontsize=12)

# FIX: Use \mathbf instead of \textbf
caption_text = (r'$\mathbf{Fig\ 1:}$ Graphical solution for 3-wave resonance.' + '\n' +
               r'Green Circle: $\Omega(\kappa) = \Omega(\kappa - \kappa_0) + \Omega(\kappa_0)$')
fig.text(0.5, 0.02, caption_text, ha='center', fontsize=11, 
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

# %% 4. SLIDER SETUP

ax_slider = fig.add_subplot(gs[1])
ax_slider.axis('off') 

slider_ax_rect = plt.axes([0.2, 0.1, 0.6, 0.03]) 
slider = Slider(
   ax=slider_ax_rect,
   label=r'Adjust $\kappa_0$: ',
   valmin=k0_min,
   valmax=k0_max,
   valinit=k0_init,
   valstep=k0_step,
   color='gray'
)

val_text = fig.text(0.85, 0.1, r'$\kappa_0 = 0.00$', fontsize=14)

# %% 5. UPDATE LOGIC

def update(val):
   k0 = slider.val
   val_text.set_text(r'$\kappa_0 = %.2f$' % k0)
   w0 = get_omega(k0, Bo)
   
   point_k0.set_data([k0], [w0])
   
   mask = k_grid >= k0
   x_local = k_grid[mask]
   k_input = x_local - k0
   
   omega_local = get_omega(k_input, Bo)
   omega_shifted = omega_local + w0
   
   line_shift.set_data(x_local, omega_shifted)
   
   y_base_segment = omega_base[mask]
   diff = y_base_segment - omega_shifted
   
   crossings = np.where(diff[:-1] * diff[1:] <= 0)[0]
   valid_crossings = [idx for idx in crossings if x_local[idx] > k0 + 0.1]
   
   if valid_crossings:
       idx = valid_crossings[0] 
       x1, x2 = x_local[idx], x_local[idx+1]
       y1, y2 = diff[idx], diff[idx+1]
       
       # FIX: Added small epsilon (1e-12) to prevent division by zero
       denom = abs(y1) + abs(y2)
       fraction = abs(y1) / (denom + 1e-12)
       
       k_int_exact = x1 + fraction * (x2 - x1)
       w_int_exact = np.interp(k_int_exact, k_grid, omega_base)
       
       point_int.set_data([k_int_exact], [w_int_exact])
       point_int.set_visible(True)
   else:
       point_int.set_visible(False)
       
   fig.canvas.draw_idle()

slider.on_changed(update)
update(k0_init)

plt.show()