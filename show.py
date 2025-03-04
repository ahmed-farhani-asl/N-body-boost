import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.animation import FuncAnimation

#G = 39.2615
dt = 0.005

# Initialize the figure and 3D plot
fig = plt.figure(figsize=(14, 7))
ax_com = fig.add_subplot(121, projection='3d')  # 3D plot with center of mass at origin
ax_energy = fig.add_subplot(122)  # Energy plot

plt.subplots_adjust(wspace=0.5)

# Labels for the 3D plot
ax_com.set_xlabel('X [AU]', fontsize=12)
ax_com.set_ylabel('Y [AU]', fontsize=12)
ax_com.set_zlabel('Z [AU]', fontsize=12)

# Energy plot labels
ax_energy.set_xlabel('Time Step')
ax_energy.set_ylabel('Energy')

# Read the CSV file
csv_file = "simulation_output.csv"
df = pd.read_csv(csv_file)

# Extract number of bodies and positions from the CSV
num_bodies = (df.shape[1] - 4) // 6  
positions = np.array([df.iloc[:, (4 + i * 6):(4 + i * 6 + 3)].values for i in range(num_bodies)])
positions = positions.transpose(1, 0, 2)  # Shape: (timesteps, bodies, 3)

all_frames = len(positions)

# Extract energy data
energy_data = {'KE': df['KE'].values, 'PE': df['PE'].values, 'Total': df['Total_Energy'].values, 'dE': np.array([(U - df['Total_Energy'].values[0])/U for U in df['Total_Energy'].values])}

# Initialize scatter plot
#colors = ["red", "blue", "green", "brown", "orange", "purple", "pink", "brown", "black", "cyan"]
scatter = ax_com.scatter(positions[0, :, 0], positions[0, :, 1], positions[0, :, 2], c='g', marker='o')

# Initialize trace lines
trace_lines = [ax_com.plot([], [], [], c='r', lw=1, alpha=0.5)[0] for i in range(num_bodies)]

speed = 20

# Update function for animation
def update(frame):

    frame *= speed
    
    if frame >= all_frames:
        return scatter, *trace_lines

    # Update scatter plot
    scatter._offsets3d = (positions[frame, :, 0], positions[frame, :, 1], positions[frame, :, 2])

    # Update trace lines
    for i in range(num_bodies):
        start = max(0, frame - 5000)
        trace_lines[i].set_data(positions[start:frame, i, 0], positions[start:frame, i, 1])
        trace_lines[i].set_3d_properties(positions[start:frame, i, 2])

    # Update dynamic axis limits
    xlimit = [positions[frame, :, 0].min() - 10, positions[frame, :, 0].max() + 10]
    ylimit = [positions[frame, :, 1].min() - 10, positions[frame, :, 1].max() + 10]
    zlimit = [positions[frame, :, 2].min() - 10, positions[frame, :, 2].max() + 10]
#    xlimit = [-40, 40]
#    ylimit = [-40, 40]
#    zlimit = [-15, 15]
    ax_com.set_xlim(xlimit)
    ax_com.set_ylim(ylimit)
    ax_com.set_zlim(zlimit)
    ax_com.set_xticks(xlimit)
    ax_com.set_yticks(ylimit)
    ax_com.set_zticks(zlimit)
    ax_com.set_aspect('equal')
    ax_com.grid(True)

    # Update energy plot
    ax_energy.clear()
    ax_energy.set_xlabel('Time [year]', fontsize=14)
    ax_energy.set_ylabel('Energy', fontsize=14)
    ax_energy.plot(np.array([i * dt for i in range(1, frame + 2)]), energy_data['KE'][:frame + 1], label='Kinetic Energy (KE)', color='r')
    ax_energy.plot(np.array([i * dt for i in range(1, frame + 2)]), energy_data['PE'][:frame + 1], label='Potential Energy (PE)', color='b')
    ax_energy.plot(np.array([i * dt for i in range(1, frame + 2)]), energy_data['Total'][:frame + 1], label='Total Energy (UE)', color='g')
    ax_energy.legend(loc='upper right')

    return scatter, *trace_lines

# Create the animation
ani = FuncAnimation(fig, update, frames=all_frames-speed-1, interval=1, blit=False)

# Show the animation
plt.show()

#ani.save('animation7.gif', fps=5, dpi=70)
#plt.savefig('final_figure.png', dpi=300)