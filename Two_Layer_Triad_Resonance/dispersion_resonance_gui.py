import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox

# =============================================================================
# 1. USER CONFIGURATION (EDIT HERE)
# =============================================================================
CONFIG = {
    # --- Initial Simulation Values ---
    'init': {
        'k0': 1.5,      # Shift Point (kappa_0)
        'R':  0.85,     # Density Ratio (rho1/rho2)
        'h':  1.0       # Depth Ratio (h1/h2)
    },
    
    # --- Initial Plot Limits ---
    'limits': {
        'xlim': 8.0,    # Plot Range X (+/- kappa)
        'ymax': 4.0,    # Y Max (Omega)
        'ymin': 0.0     # Y Min (Omega)
    },

    # --- Slider Ranges [Min, Max, Step] ---
    'ranges': {
        'k0':   [0.1,   5.0,   0.01],
        'R':    [0.01,  0.99,  0.001],
        'h':    [0.1,   5.0,   0.01],
        'xlim': [2.0,   20.0,  0.1],
        'ymax': [0.1,   15.0,  0.1],
        'ymin': [-15.0, 5.0,   0.1]
    },

    # --- Visual Styles ---
    'style': {
        'orig_surf':  {'color': 'b', 'linestyle': '-',  'linewidth': 3},
        'orig_int':   {'color': 'r', 'linestyle': '--', 'linewidth': 3},
        'shift_surf': {'color': 'b', 'linestyle': '-',  'linewidth': 1.0},
        'shift_int':  {'color': 'r', 'linestyle': ':',  'linewidth': 1.5}
    }
}

# =============================================================================
# 2. MATH KERNEL
# =============================================================================
def solve_dispersion(k_in, R, h):
    k = np.array(k_in)
    k[np.abs(k) < 1e-9] = 1e-9
    
    th_k = np.tanh(k)
    th_kh = np.tanh(k * h)
    coth_k = 1.0 / th_k
    coth_kh = 1.0 / th_kh
    
    A = R + (coth_k * coth_kh)
    B = -k * (coth_k + coth_kh)
    C = (k**2) * (1 - R)
    
    discriminant = np.sqrt(B**2 - 4 * A * C)
    
    Omega2_1 = (-B + discriminant) / (2 * A)
    Omega2_2 = (-B - discriminant) / (2 * A)
    
    Omega2_1[Omega2_1 < 0] = 0
    Omega2_2[Omega2_2 < 0] = 0
    
    Om_Surf = np.sqrt(Omega2_1)
    Om_Int  = np.sqrt(Omega2_2)
    
    # Fix k=0
    mask_zero = np.abs(k_in) < 1e-8
    Om_Surf[mask_zero] = 0
    Om_Int[mask_zero]  = 0
    
    return Om_Surf, Om_Int

def find_intersections(x_vec, y1_vec, y2_vec):
    diff = y1_vec - y2_vec
    sign_changes = np.where(np.diff(np.sign(diff)))[0]
    
    x_int = []
    y_int = []
    
    for idx in sign_changes:
        x_a, x_b = x_vec[idx], x_vec[idx+1]
        d_a, d_b = diff[idx], diff[idx+1]
        
        if d_a == d_b: continue 
            
        frac = abs(d_a) / (abs(d_a) + abs(d_b))
        x_val = x_a + frac * (x_b - x_a)
        
        val_a, val_b = y1_vec[idx], y1_vec[idx+1]
        y_val = val_a + frac * (val_b - val_a)
        
        x_int.append(x_val)
        y_int.append(y_val)
        
    return x_int, y_int

# =============================================================================
# 3. GUI SETUP
# =============================================================================
class ControlRow:
    def __init__(self, fig, y_pos, label_latex, key, update_func):
        # Unpack config for this specific key
        val_min, val_max, step = CONFIG['ranges'][key]
        
        # Determine initial value based on key type
        if key in CONFIG['init']:
            val_init = CONFIG['init'][key]
        else:
            val_init = CONFIG['limits'][key]

        self.step = step
        self.val_min = val_min
        self.val_max = val_max
        self.update_func = update_func

        # Label
        fig.text(0.02, y_pos + 0.055, label_latex, fontsize=11, fontweight='bold', ha='left')

        # Slider
        ax_slider = plt.axes([0.02, y_pos, 0.20, 0.03])
        # We allow default formatting initially to avoid TypeError
        self.slider = Slider(ax_slider, '', val_min, val_max, valinit=val_init, valstep=step)
        
        # --- FIX: Hide the number outside the slider manually ---
        self.slider.valtext.set_visible(False)
        
        # TextBox
        ax_box = plt.axes([0.23, y_pos, 0.08, 0.04])
        self.textbox = TextBox(ax_box, '', initial=str(val_init))

        self.slider.on_changed(self.on_slider_change)
        self.textbox.on_submit(self.on_text_submit)

    def on_slider_change(self, val):
        self.textbox.eventson = False 
        self.textbox.set_val(f"{val:.3f}")
        self.textbox.eventson = True
        self.update_func(None)

    def on_text_submit(self, text):
        try:
            val = float(text)
            if val < self.val_min: val = self.val_min
            if val > self.val_max: val = self.val_max
            
            self.slider.eventson = False
            self.slider.set_val(val)
            self.slider.eventson = True
            self.update_func(None)
        except ValueError:
            self.textbox.set_val(str(self.slider.val))

    @property
    def val(self):
        return self.slider.val

# Create Figure
fig = plt.figure(figsize=(14, 9))
plt.subplots_adjust(left=0.40, bottom=0.10, right=0.95, top=0.95)

# Axis Setup
ax = fig.add_subplot(111)
ax.set_xlabel(r'$\kappa$', fontsize=20)
ax.set_ylabel(r'$\Omega$', fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.grid(True)

# Prepare Line Objects using Config Styles
lines = {}
lines['orig_s1'], = ax.plot([], [], **CONFIG['style']['orig_surf'])
lines['orig_s2'], = ax.plot([], [], **CONFIG['style']['orig_surf'])
lines['orig_i1'], = ax.plot([], [], **CONFIG['style']['orig_int'])
lines['orig_i2'], = ax.plot([], [], **CONFIG['style']['orig_int'])

lines['shift_s1'], = ax.plot([], [], **CONFIG['style']['shift_surf'])
lines['shift_s2'], = ax.plot([], [], **CONFIG['style']['shift_surf'])
lines['shift_i1'], = ax.plot([], [], **CONFIG['style']['shift_int'])
lines['shift_i2'], = ax.plot([], [], **CONFIG['style']['shift_int'])

# --- MODIFIED MARKERS ---
# 1. Origin Point (Fixed Black Circle)
pt_origin, = ax.plot([0], [0], 'ko', markersize=11, zorder=11)

# 2. Shift Center Point (Black Circle)
pt_center, = ax.plot([], [], 'ko', markersize=11, zorder=11)

# 3. Other Intersections (Yellow Circles, NO black edge)
# Changed markeredgecolor='k' to markeredgecolor='none'
# Use 'o' for circle, pure Hex Yellow '#FFFF00', and size 12 (30% bigger)
pt_intersect, = ax.plot([], [], 'o', color='#FFFF00', markersize=12, markeredgecolor='k', zorder=10)

# Legend Info
info_text = (
    "Thick Blue: Orig. Surface\n"
    "Thick Red:  Orig. Internal\n"
    "Thin Blue:  Shifted Surface\n"
    "Thin Red:   Shifted Internal\n"
    "Black O:    Origin & Shift Point\n"
    "Yellow O:   Other Intersections"
)
fig.text(0.02, 0.15, info_text, fontsize=10, style='italic', va='top', bbox=dict(facecolor='white', alpha=0.5))
fig.text(0.02, 0.96, "Control Panel", fontsize=14, fontweight='bold')

# =============================================================================
# 4. UPDATE LOGIC
# =============================================================================
def update(val):
    xl = ctrl_xlim.val
    yl_max = ctrl_ymax.val
    yl_min = ctrl_ymin.val
    k0 = ctrl_k0.val
    R = ctrl_R.val
    h = ctrl_h.val

    if yl_min >= yl_max: yl_min = yl_max - 0.1

    calc_limit = xl * 1.1
    k_vec = np.linspace(-calc_limit, calc_limit, 1000)

    # 1. Original
    Om_S_orig, Om_I_orig = solve_dispersion(k_vec, R, h)
    Orig_Y = [Om_S_orig, -Om_S_orig, Om_I_orig, -Om_I_orig]

    # 2. Shift Point
    Om_S_0, _ = solve_dispersion([k0], R, h)
    Omega0 = Om_S_0[0]

    # 3. Shifted
    k_shift = k_vec - k0
    Om_S_shift, Om_I_shift = solve_dispersion(k_shift, R, h)
    Shift_Y = [
        Omega0 + Om_S_shift, Omega0 - Om_S_shift,
        Omega0 + Om_I_shift, Omega0 - Om_I_shift
    ]

    # Update Lines
    lines['orig_s1'].set_data(k_vec, Orig_Y[0])
    lines['orig_s2'].set_data(k_vec, Orig_Y[1])
    lines['orig_i1'].set_data(k_vec, Orig_Y[2])
    lines['orig_i2'].set_data(k_vec, Orig_Y[3])

    lines['shift_s1'].set_data(k_vec, Shift_Y[0])
    lines['shift_s2'].set_data(k_vec, Shift_Y[1])
    lines['shift_i1'].set_data(k_vec, Shift_Y[2])
    lines['shift_i2'].set_data(k_vec, Shift_Y[3])

    # Update Shift Center Point (Black)
    pt_center.set_data([k0], [Omega0])

    # Find Intersections
    ix_vals = []
    iy_vals = []

    for y1 in Orig_Y:
        for y2 in Shift_Y:
            tx, ty = find_intersections(k_vec, y1, y2)
            for x_i, y_i in zip(tx, ty):
                if yl_min <= y_i <= yl_max:
                    # Filter out the origin and the shift point from the yellow list
                    tol = 1e-5
                    is_origin = (abs(x_i) < tol and abs(y_i) < tol)
                    is_center = (abs(x_i - k0) < tol and abs(y_i - Omega0) < tol)
                    
                    if not is_origin and not is_center:
                        ix_vals.append(x_i)
                        iy_vals.append(y_i)

    # Update Yellow Intersections
    pt_intersect.set_data(ix_vals, iy_vals)
    
    ax.set_xlim(-xl, xl)
    ax.set_ylim(yl_min, yl_max)
    fig.canvas.draw_idle()

# =============================================================================
# 5. INITIALIZE CONTROLS
# =============================================================================
curr_y = 0.85
def next_y():
    global curr_y
    y = curr_y
    curr_y -= 0.10
    return y

# Generate Controls (Order matters for layout)
ctrl_xlim = ControlRow(fig, next_y(), r'Plot Range $X (\pm\kappa)$', 'xlim', update)
ctrl_ymax = ControlRow(fig, next_y(), r'Y Max ($\Omega$)', 'ymax', update)
ctrl_ymin = ControlRow(fig, next_y(), r'Y Min ($\Omega$)', 'ymin', update)

# Divider
line_ax = plt.axes([0.02, curr_y + 0.05, 0.30, 0.002])
line_ax.set_xticks([]); line_ax.set_yticks([]); line_ax.set_facecolor('black')

ctrl_k0 = ControlRow(fig, next_y(), r'Shift Point $(\kappa_0)$', 'k0', update)
ctrl_R  = ControlRow(fig, next_y(), r'Density Ratio $(R)$', 'R', update)
ctrl_h  = ControlRow(fig, next_y(), r'Depth Ratio $(h)$', 'h', update)

update(None)
plt.show()