import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox

# Define Constants
h_val = 1.0
mu_val = 1.0

# Initialize Figure
fig, ax = plt.subplots(figsize=(9, 6))
plt.subplots_adjust(left=0.1, bottom=0.35)

ax.set_title('Velocity Profile $u(y)$')
ax.set_xlabel('Velocity $u$ (m/s)')
ax.set_ylabel('Channel Height $y$ (m)')
ax.grid(True)
ax.axhline(0, color='k', lw=3)
ax.axhline(h_val, color='k', lw=3)
ax.set_xlim(-15, 15)
ax.set_ylim(0, h_val)

# Define grid and extraction indices for arrows
y = np.linspace(0, h_val, 100)
idx = np.arange(0, len(y), 6) # Spacing to prevent vertical overlap
y_q = y[idx]

# Initialize Plot Objects
line, = ax.plot([], [], 'b-', lw=2)

# Quiver settings optimized for small, filled, exact-length arrowheads
quiv = ax.quiver(np.zeros_like(y_q), y_q, np.zeros_like(y_q), np.zeros_like(y_q),
                color='r', angles='xy', scale_units='xy', scale=1,
                width=0.003, headwidth=4, headlength=5, headaxislength=4)

# Create UI Axes
ax_U_sld = plt.axes([0.25, 0.2, 0.45, 0.03])
ax_U_txt = plt.axes([0.75, 0.2, 0.1, 0.03])
ax_P_sld = plt.axes([0.25, 0.1, 0.45, 0.03])
ax_P_txt = plt.axes([0.75, 0.1, 0.1, 0.03])

# Create UI Widgets
sld_U = Slider(ax_U_sld, 'Top Velocity ($U$)', -10.0, 10.0, valinit=0.0)
txt_U = TextBox(ax_U_txt, '', initial='0.00')
sld_P = Slider(ax_P_sld, 'Pressure Grad ($dP/dx$)', -20.0, 20.0, valinit=0.0)
txt_P = TextBox(ax_P_txt, '', initial='0.00')

# State flag to prevent recursion between text and slider updates
updating = False

def update_plot():
   U = sld_U.val
   dPdx = sld_P.val
   
   # Calculate profile
   term1 = (y / h_val) * U
   term2 = (h_val**2 / (2 * mu_val)) * (-dPdx) * (y / h_val) * (1 - (y / h_val))
   u = term1 + term2
   
   line.set_data(u, y)
   
   # Update arrows and hide near-zero vectors
   u_q = u[idx]
   mask = np.abs(u_q) > 0.05
   u_q_masked = np.where(mask, u_q, 0)
   quiv.set_UVC(u_q_masked, np.zeros_like(u_q_masked))
   
   fig.canvas.draw_idle()

def sync_U_from_sld(val):
   global updating
   if not updating:
       updating = True
       txt_U.set_val(f"{val:.2f}")
       update_plot()
       updating = False

def sync_P_from_sld(val):
   global updating
   if not updating:
       updating = True
       txt_P.set_val(f"{val:.2f}")
       update_plot()
       updating = False

def sync_U_from_txt(text):
   global updating
   if not updating:
       updating = True
       try:
           val = max(min(float(text), 10.0), -10.0)
           sld_U.set_val(val)
           txt_U.set_val(f"{val:.2f}")
       except ValueError:
           txt_U.set_val(f"{sld_U.val:.2f}")
       update_plot()
       updating = False

def sync_P_from_txt(text):
   global updating
   if not updating:
       updating = True
       try:
           val = max(min(float(text), 20.0), -20.0)
           sld_P.set_val(val)
           txt_P.set_val(f"{val:.2f}")
       except ValueError:
           txt_P.set_val(f"{sld_P.val:.2f}")
       update_plot()
       updating = False

# Connect Listeners
sld_U.on_changed(sync_U_from_sld)
sld_P.on_changed(sync_P_from_sld)
txt_U.on_submit(sync_U_from_txt)
txt_P.on_submit(sync_P_from_txt)

# Initial Draw
update_plot()
plt.show()
