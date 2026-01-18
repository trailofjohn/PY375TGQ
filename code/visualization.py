
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider, Button
import numpy as np
from simulation_core import AU, DAY

def run_visualization(hist, title):
    print("Generating Animation...")
    plt.style.use('dark_background')
    
    fig = plt.figure(figsize=(18, 10))
    fig.suptitle(f"Spacetime Estimate: {title}", color='white', fontsize=16, fontweight='bold')
    
    # Layout:
    # [ 3D Trajectory (Large) ] [ Relative Error (Top Right) ]
    # [                       ] [ Clock Bias (Mid Right)     ]
    # [                       ] [ Shapiro Delay (Bot Right)  ]
    
    gs = fig.add_gridspec(3, 3, width_ratios=[2, 1, 1], hspace=0.5)
    ax3d = fig.add_subplot(gs[:, 0], projection='3d')
    ax_err = fig.add_subplot(gs[0, 1:]) # Top Right spanning 2 cols
    ax_clk = fig.add_subplot(gs[1, 1:])
    ax_shap = fig.add_subplot(gs[2, 1:])
    
    # Widgets Area
    # Slider Axis
    ax_slider = plt.axes([0.25, 0.02, 0.35, 0.03], facecolor='black')
    # Speed: 0.1x to 10x. Logarithmic feel by using low linear range.
    speed_slider = Slider(ax_slider, 'Speed', 0.1, 10.0, valinit=1.0, valstep=0.1, color='cyan')
    
    # Play/Pause Button
    ax_play = plt.axes([0.1, 0.02, 0.1, 0.04])
    btn_play = Button(ax_play, 'Pause', color='black', hovercolor='dimgray')
    btn_play.label.set_color('white')
    
    anim_running = True
    def toggle_play(event):
        nonlocal anim_running
        if anim_running:
            ani.event_source.stop()
            btn_play.label.set_text("Play")
            anim_running = False
        else:
            ani.event_source.start()
            btn_play.label.set_text("Pause")
            anim_running = True
    
    btn_play.on_clicked(toggle_play)
    
    # Time Display
    time_text = fig.text(0.02, 0.95, '', color='white', fontsize=14, fontfamily='monospace')
    
    # --- 3D Plot ---
    ax3d.set_xlabel('X (AU)')
    ax3d.set_ylabel('Y (AU)')
    ax3d.set_zlabel('Z (AU)')
    ax3d.xaxis.set_pane_color((0,0,0,0)); ax3d.yaxis.set_pane_color((0,0,0,0)); ax3d.zaxis.set_pane_color((0,0,0,0))
    ax3d.grid(False)
    
    # Starfield
    stars = np.random.uniform(-5, 5, (150, 3))
    ax3d.scatter(stars[:,0], stars[:,1], stars[:,2], c='white', s=0.3, alpha=0.4)
    ax3d.scatter([0], [0], [0], color='yellow', s=300, label='Sun', edgecolors='orange', alpha=0.9)
    
    line_true, = ax3d.plot([], [], [], 'g-', label='True', linewidth=1.0, alpha=0.6)
    line_est, = ax3d.plot([], [], [], 'c:', label='EKF', linewidth=1.5)
    point_est, = ax3d.plot([], [], [], 'co', markersize=4, markeredgecolor='white')
    
    # Ellipsoid (3 Sigma) - Placeholder, will update in loop
    ellipsoid_wire = ax3d.plot_wireframe(np.array([[]]), np.array([[]]), np.array([[]]), color='cyan', alpha=0.2, linewidth=0.5)

    ax3d.legend(loc='upper left', facecolor='black', edgecolor='white', fontsize='small')
    
    limit = np.max(np.abs(hist['true_pos'])) / AU * 1.1
    ax3d.set_xlim(-limit, limit); ax3d.set_ylim(-limit, limit); ax3d.set_zlim(-limit/2, limit/2)

    # --- Relative Error Plot (The "Zoom") ---
    # Shows Position Error Magnitude over time
    ax_err.set_title("Position Error Magnitude (m) - Log Scale")
    ax_err.set_xlabel("Time (Days)")
    ax_err.set_yscale("log")
    ax_err.grid(True, alpha=0.2, which='both')
    
    # Calculate errors
    pos_err = np.linalg.norm(hist['est_pos'] - hist['true_pos'], axis=1)
    dr_err = np.linalg.norm(hist['dead_reckoning'] - hist['true_pos'], axis=1)
    
    line_err_ekf, = ax_err.plot([], [], 'c-', label='EKF Error', lw=1.5)
    line_err_dr, = ax_err.plot([], [], 'r--', label='Dead Reckoning', lw=1, alpha=0.7)
    ax_err.legend(loc='upper right', fontsize='x-small')
    
    # --- Clock Bias ---
    ax_clk.set_title("Clock Bias Residual (m)")
    ax_clk.grid(True, alpha=0.2)
    line_clk, = ax_clk.plot([], [], 'c-', lw=1)
    
    # --- Shapiro ---
    ax_shap.set_title("Shapiro Delay (m)")
    ax_shap.set_xlabel("Time (days)")
    ax_shap.grid(True, alpha=0.2)
    line_shap, = ax_shap.plot([], [], 'violet', lw=1.5)
    
    days = hist['t'] / DAY
    
    # Pre-calculate Ellipsoid spheres for efficiency if possible, or usually just re-calc
    u = np.linspace(0, 2 * np.pi, 10)
    v = np.linspace(0, np.pi, 10)
    x_sphere = np.outer(np.cos(u), np.sin(v))
    y_sphere = np.outer(np.sin(u), np.sin(v))
    z_sphere = np.outer(np.ones(np.size(u)), np.cos(v))

    # We need to track our own "virtual frame" because the slider changes how fast we move through history
    current_idx = 0.0
    
    def update(frame):
        nonlocal current_idx
        
        # Read scalar speed
        speed = speed_slider.val 
        
        # Increment virtual index
        current_idx += speed
        
        # Loop functionality or Stop at end? Let's loop for demo or stop.
        if current_idx >= len(days): 
            current_idx = 0 # Loop
        
        idx = int(current_idx)
        if idx == 0: idx = 1
        
        # Update Time Display
        time_text.set_text(f"Mission Elapsed Time: {days[idx]:.2f} Days")

        # 3D Trajectory
        # Show full trail
        line_true.set_data(hist['true_pos'][:idx, 0]/AU, hist['true_pos'][:idx, 1]/AU)
        line_true.set_3d_properties(hist['true_pos'][:idx, 2]/AU)
        
        line_est.set_data(hist['est_pos'][:idx, 0]/AU, hist['est_pos'][:idx, 1]/AU)
        line_est.set_3d_properties(hist['est_pos'][:idx, 2]/AU)
        
        # Current Position
        cx, cy, cz = hist['est_pos'][idx, 0]/AU, hist['est_pos'][idx, 1]/AU, hist['est_pos'][idx, 2]/AU
        point_est.set_data([cx], [cy])
        point_est.set_3d_properties([cz])
        
        # Uncertainty Ellipsoid (3-Sigma)
        # P matrix trace was saved, but ideally we need full P or at least diagonals.
        # Approximation: Using Trace/3 as radius for a sphere since we didn't save full P history
        # (Modifying core to save diagonals is expensive for history memory maybe? Let's assume sphere for demo)
        if 'P_trace' in hist:
            sigma_r = np.sqrt(hist['P_trace'][idx] / 3.0) * 3 # 3-sigma radius in meters
            sigma_r_au = sigma_r / AU
            
            # Scale sphere
            # To make it visible at solar system scale we might need to artifically scale it 
            # OR only show it when it's huge. 
            # Scaling it by a factor for visibility (e.g. 1000x) and noting it?
            # User wants to SEE it.
            visual_scale = 1000.0 # Multiplier for visibility
            
            nonlocal ellipsoid_wire
            if ellipsoid_wire: ellipsoid_wire.remove()
            ellipsoid_wire = ax3d.plot_wireframe(
                 x_sphere * sigma_r_au * visual_scale + cx,
                 y_sphere * sigma_r_au * visual_scale + cy,
                 z_sphere * sigma_r_au * visual_scale + cz,
                 color='cyan', alpha=0.3, rstride=1, cstride=1, linewidth=0.5
            )
        
        # Error Plots
        line_err_ekf.set_data(days[:idx], pos_err[:idx])
        line_err_dr.set_data(days[:idx], dr_err[:idx])
        ax_err.set_xlim(0, max(days[:idx]) + 1)
        ax_err.set_ylim(min(pos_err[pos_err>0]) * 0.5, max(dr_err[:idx]) * 1.5) # Log scale handling
        
        # Clock
        line_clk.set_data(days[:idx], hist['clock_err'][:idx])
        ax_clk.set_xlim(0, max(days[:idx]) + 1)
        current_clk_errs = hist['clock_err'][:idx]
        if len(current_clk_errs) > 0:
            mx = np.max(np.abs(current_clk_errs))
            ax_clk.set_ylim(-mx*1.2, mx*1.2)
            
        # Shapiro
        line_shap.set_data(days[:idx], hist['shapiro'][:idx])
        ax_shap.set_xlim(0, max(days[:idx]) + 1)
        ax_shap.set_ylim(0, max(hist['shapiro']) * 1.1)
        
        return line_true, line_est, point_est

    # Frames calculation
    # We want valid interaction, so we use a generator or large number of frames
    ani = animation.FuncAnimation(fig, update, frames=200, interval=20, blit=False, cache_frame_data=False)
    
    plt.show()
