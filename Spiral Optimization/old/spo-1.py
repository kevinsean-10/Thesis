import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Define the spiral-shaped function to optimize
def spiral_function(x, y):
    return (x**2 + y**2) + np.sin(np.sqrt(x**2 + y**2)) * 0.1

# Parameters for the optimization
num_particles = 1
num_iterations = 300
learning_rate = 0.1

# Initialize particles with random positions
particles = np.random.uniform(-5, 5, size=(num_particles, 2))
velocities = np.zeros_like(particles)

# Create a figure for animation
fig, ax = plt.subplots()
ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)
scatter = ax.scatter([], [], c='red', marker='o')

# Update function for animation
def update(frame):
    global particles, velocities
    for i in range(num_particles):
        gradient_x = 2 * particles[i, 0] + 0.1 * np.cos(np.sqrt(particles[i, 0]**2 + particles[i, 1]**2)) * particles[i, 0] / np.sqrt(particles[i, 0]**2 + particles[i, 1]**2)
        gradient_y = 2 * particles[i, 1] + 0.1 * np.cos(np.sqrt(particles[i, 0]**2 + particles[i, 1]**2)) * particles[i, 1] / np.sqrt(particles[i, 0]**2 + particles[i, 1]**2)
        velocities[i, 0] = learning_rate * gradient_x
        velocities[i, 1] = learning_rate * gradient_y
        particles[i] -= velocities[i]
    
    scatter.set_offsets(particles)
    return scatter,

# Create the animation
ani = FuncAnimation(fig, update, frames=num_iterations, blit=True)

# Show the animation
plt.show()
