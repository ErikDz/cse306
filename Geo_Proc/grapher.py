
import matplotlib.pyplot as plt

# Read the time taken to compute each frame from the file
with open('build/time_dt=0.002.txt', 'r') as file:
    times = [float(line.strip()) for line in file]

# Plot the time taken for each frame
plt.plot(range(1, len(times) + 1), times)
plt.xlabel('Frame')
plt.ylabel('Time (seconds)')
plt.title('Time taken to compute each frame')

# we save the graph
plt.savefig('time_taken.png')